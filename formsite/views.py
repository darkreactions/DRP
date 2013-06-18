from django.http import HttpResponse, Http404, HttpResponseRedirect
from django.template.loader import get_template
from django.template import Context
from django.shortcuts import render
from django.utils import simplejson
from models import *
from django import forms
from django.views.decorators.csrf import csrf_exempt
from validation import *
import csv
import string
from django.contrib.auth.forms import UserCreationForm
from django.contrib import auth
import datetime

#Setup Variables
current_page = 1 #1-based index.
data_per_page = 100

#Helper Function that returns all data of a specific lab group.
def get_saved_data(lab_group):
	return Data.objects.filter(lab_group=lab_group)

def data_view(request, num = current_page): #If no data is entered, stay on the current page.
	#Call the global variables.
	global data_per_page
	global current_page
	
	#Make sure the user is logged in before they see the form.
	if request.user.is_authenticated():
		#Get the Lab_Group data if the user is logged in.
		u = request.user
		saved_data = get_saved_data(u.get_profile().lab_group)
		
		# If the form has been submitted.
		if request.method == 'POST': 
			#Bind the user's data and verify that it is legit.
			data_form = DataEntryForm(user=u, data=request.POST)
			if data_form.is_valid():
				#If all data is valid, save the entry (and the submitter)
				data_form.save()
		else:
			#Submit a blank form if one was not just submitted.
			data_form = DataEntryForm()
	else:
		#If the user is not logged in. ###
		saved_data = []
		data_form = DataEntryForm()
	
	#Calculate the number of possible pages.
	total_pages = int(1 + (len(saved_data)-1)/data_per_page) ###
	if (total_pages < 1): total_pages = 1
	
	###if (request.method == 'POST'): #Jump to the new data if relevant.
		###if form.is_valid(): num = total_pages 

	try:
		num = int(num) #Convert the unicode page number to an integer for comparison.
		assert 0 < num <= total_pages
	except:
		num = total_pages
	
	#Update current_page if it was changed.
	current_page = num
	
	#Only send the data on the requested page. Note that Lab_Group.saved_data is a 0-based index.
	start_index = (current_page-1)*data_per_page 
	end_index = (current_page)*data_per_page
	index_range = range(start_index+1, end_index+1) #1-based index for visuals

	data_package = zip(saved_data[start_index:end_index], index_range)
	pageLinks = get_pageLinks(current_page, total_pages, data_per_page)
	page_package = {
		"current_page":current_page,
		"pageLinks":pageLinks,
		"total_pages":total_pages,
		"data_per_page":data_per_page,
		}
		
	
	return render(request, 'database_global.html', {
		"form": data_form,
		"data_on_page": data_package, #Includes data indexes
		"page_package": page_package, #Includes data indexes
		"total_data_size": len(saved_data),
		#"user_package": user_package
	})

#Helper function that returns a list of page numbers/ellipses.
#	Note: All vars relate to "_page" (eg, "current_page")
def get_pageLinks(current, total, data_per):
	#Setup Variables
	current_radius = 4 #Number of links to display "around" current.
	
	#Always display the first page.
	pageLinks = {1} #Use a set to remove any duplicates. 
	
	if total > 1:
		for i in range(current-current_radius, current+current_radius+1):
			if (i > 1) and (i < total): pageLinks.add(i)
			
		#Always display the last page if applicable.
		pageLinks.add(total)
		
	#Convert pageLinks to an ordered list.
	pageLinks = list(pageLinks)
	pageLinks.sort()
	
	if len(pageLinks) >= 2:
		i=0
		while i < len(pageLinks):
			try:
				if pageLinks[i+1]-pageLinks[i] > 1:
					pageLinks.insert(i+1,"...")
					i+=1 #Extra addition to account for new "..." element.
			except:
				pass #At the last element of the list or an error, just skip.
			i += 1
	return pageLinks

#Helper function that returns a related data entry field ("reactant 1 name" --> "reactant_1")
def get_related_field(heading):
	#Strip all punctuation, capitalization and spacing from the header. 
	#Note, translate() is a super fast version of replace())
	stripped_heading = heading.translate(None, string.punctuation)
	stripped_heading = stripped_heading.translate(None, " ").lower()
	stripped_heading = stripped_heading[:20] #Limit the checked heading (saves time if super long).
	if ("reacta" in stripped_heading or "mass" in stripped_heading
		or "unit" in stripped_heading or "vol" in stripped_heading
		or "amou" in stripped_heading or "name" in stripped_heading
		or "qua" in stripped_heading):
		if ("mass" in stripped_heading or "quantity" in stripped_heading
			or "vol" in stripped_heading or "amount" in stripped_heading):
			related_field = "quantity_"
		elif "unit" in stripped_heading:
			related_field = "unit_"
		else:
			related_field = "reactant_"
		for i in range(5): #ie, 1-5
			if str(i+1) in stripped_heading: 
				related_field += str(i+1)
				break; #Only add 1 number to the data form.
	elif "temp" in stripped_heading:
		related_field = "temp"
	elif "time" in stripped_heading:
		related_field = "time" 
	elif "cool" in stripped_heading or "slow" in stripped_heading:
		related_field = "slow_cool"
	elif "pur" in stripped_heading:
		related_field = "purity" 
	elif "leak" in stripped_heading or "error" in stripped_heading:
		related_field = "leak" 
	elif ("ref" in stripped_heading or "cont" in stripped_heading
		or "num" in stripped_heading):
		related_field = "ref" 
	elif "out" in stripped_heading or "res" in stripped_heading:
		related_field = "outcome" 
	elif ("note" in stripped_heading or "other" in stripped_heading 
		or "info" in stripped_heading):
		related_field = "notes" 
	elif "ph" in stripped_heading:
		related_field = "pH" 
	else: #ie, related_field is unchanged.
		raise Exception("TRANSLATION ERROR: No relation found for {}.".format(heading))###Possible Raise?
	return related_field
		
def upload_data(request):
	true_fields = get_data_field_names()
	not_required = { ###Auto-generate?
			"reactant_3", "quantity_3", "unit_3", 
			"reactant_4", "quantity_4", "unit_4",
			"reactant_5", "quantity_5", "unit_5",
			"notes"
		}
			
	u = request.user
	added_quantity=0###
	error_quantity=0###
	
	if request.method=="POST" and request.FILES and u.is_authenticated():
		uploaded_file = request.FILES["file"]
		if len(uploaded_file) == 0:
			return HttpResponse("No data selected!")
		try:
			blacklist = {"x", "-1", -1, "z", "?", "", " "} #Implies absence of data. ###
			unknown_label = "UNKNOWN" #Label that blacklist values receive.
			
			#Attempt to validate the headings of the uploaded doc.
			headings_valid = False
			validation_attempt = 0
			#Separate data into groups of fields -- then separate fields.
			for data_group in csv.reader(uploaded_file, delimiter=","):
				#The first data_group should be a series of headings.
				if validation_attempt > 5: raise Exception("Unable to validate headings.")
				try:
					if not headings_valid: #Remember which column has which heading.
						#Translate the user's field name into a usable field name.
						user_fields = [get_related_field(field) for field in data_group]
						
						#If unit columns were not supplied, add them after each mass.
						set_user_fields = set(user_fields)
						auto_added_fields = set()
						for i in range(1,6): #Since only 5 reactants supported...###
							if not "unit_{}".format(i) in set_user_fields:
								user_fields.insert(
									user_fields.index("quantity_{}".format(i))+1,
									"unit_{}".format(i)
									)
									###Check if mass is last field (fence-post?).
									###Should be fine, but check anyway. ; )
								auto_added_fields.add("unit_{}".format(i))
								
						#Assert that there are no duplicates in the list
						#	and that all fields are valid.
						assert(len(user_fields)==len(true_fields))
						for field in true_fields:
							assert(field in user_fields)
							
						headings_valid = True
						continue
				except Exception as e:
					print e
					validation_attempt += 1
					continue
					
				#All other rows should have data corresponding to the headings.
				try:
					#Create new object that will receive the uploaded data.
					entry_fields = {}
					
					#Access the data from the uploaded file.
					i = 0 # data_group index (gets a datum)
					j = 0 # user_fields index. (gets a field)
					while (i < len(data_group) or j < len(user_fields)):
						#Required since data and fields may be disjunct from missing units.
						datum = data_group[i]
						field = user_fields[j]
					
						if field in auto_added_fields: 
							#Skip the field if it was auto-added because
							#	no data is present in the generated column.
							j += 1
							continue
						try:
							#If the datum isn't helpful, don't remember it.###
							if datum in blacklist: 
								if field in not_required:
									datum = "" #Mark the field as a blank (make unknown data uniform).
								else:
									datum = unknown_label ###Take this value or no?
							else:
								#If the field is a quantity, check for units.
								if field[:-2]=="quantity":
									#Remove punctuation and whitespace if necessary.
									datum = str(datum).lower()
									datum = datum.translate(None, "\n?/,!@#$%^&*-+=_\\|") #Remove gross stuff.
									stripped_datum = ""
									unit = "" #Gather a unit from quantity if present.
									for element in datum:
										if not element in "1234567890.":
											datum = float(stripped_datum)
											
											#If another element is reached, assume it is a unit. 
											#	(if it isn't valid, raise an exception)
											if element == "g": unit = "g"
											if element == " ": unit = "g" #Mark data before a space as grams.
											elif element == "m": unit = "ml"
											elif element == "d": unit = "d"
											else: raise Exception("Unknown unit in quantity.")	
											break #Ignore anything beyond the unit.
										else:
											stripped_datum += element
											
									#If the unit was auto-added_quantity, add the unit to the correct field.
									corresponding_unit = "unit_{}".format(field[-1])
									if corresponding_unit in auto_added_fields:
										#If no unit was gathered, default to grams.
										if unit=="": unit = "g"
										#Apply the unit to the data entry.
										entry_fields[corresponding_unit] = unit
								try:
									#Attempt to validate the data.
									if datum != unknown_label:
										assert(quick_validation(field, datum))
								except:
									raise Exception("Data did not pass validation!")
							entry_fields[field] = datum
							#Continue to iterate through the data_group.
							i+=1
							j+=1
						except:
							raise Exception("Entry could not be added!")
					
					#Add the new entry to the database. ###SLOWWwwwww...
					create_data_entry(u, **entry_fields)
					added_quantity += 1###
				except Exception as e:
					print "ERROR:",e,"\n",data_group,"at",field,"({})".format(datum)
					error_quantity +=1###
					
		except Exception as e:
			print "ERROR:", e
	
	#Bulk Upload?
	
	
	message = "{} entries added.".format(added_quantity)
	message += "\n{} entries failed validation.".format(error_quantity)
	return HttpResponse(message)

def download_CSV(request):
	u = request.user
	if u.is_authenticated():
		#Generate a file name.
		date = datetime.datetime.now()
		file_name = "{:2}_{:0>2}_{:0>2}_{}".format(u.get_profile().lab_group.lab_title,###
			date.day, date.month, date.year)
		
		CSV_file = HttpResponse(content_type="text/csv")
		CSV_file["Content-Disposition"] = "attachment; filename={}.csv".format(file_name)
		
		#Django HttpResponse objects can be handleded like files.
		writer = csv.writer(CSV_file)
		
		#Write the verbose headers to the CSV_file
		verbose_headers = get_data_field_names(True)
		writer.writerow(verbose_headers)
		
		
		#Write the actual entries to the CSV_file if the user is authenticated.
		headers = get_data_field_names()
		errors_total = 0###TEST VARS
		#Get the Lab_Group data to allow direct manipulation.
		saved_data = get_saved_data(u.get_profile().lab_group)	
		
		for entry in saved_data:
			row = []
			try:
				for field in headers:
					row += [eval("entry.{}".format(field))]
				writer.writerow(row)
			except:
				errors_total += 1
				print "ERROR AT:\n", entry###
		print "Total errors: {}".format(errors_total)
	return CSV_file #ie, return HttpResponse(content_type="text/csv")
	
def data_transmit(request, num):
	global data_per_page
	try:
		if request.method == "GET" and request.user.is_authenticated():
			#Get the Lab_Group data if the user is logged in.
			u = request.user
			saved_data = get_saved_data(u.get_profile().lab_group)
			
			total_data_size = len(saved_data)
			total_pages = int(1 + (total_data_size-1)/data_per_page) ###
			if (total_pages < 1): total_pages = 1
			try:
				#current_page == num 
				num = int(num) #Convert the unicode page number to an integer for comparison.
				assert 0 < num <= total_pages
			except:
				print "Can't go there!"
				num = total_pages
			
			#Update current_page if it was changed.
			current_page = num
			
			#Only send the data on the requested page. Note that Lab_Group.saved_data is a 0-based index.
			start_index = (current_page-1)*data_per_page 
			end_index = (current_page)*data_per_page
			index_range = range(start_index+1, end_index+1) #1-based index for visuals

			data_package = zip(saved_data[start_index:end_index], index_range)
			pageLinks = get_pageLinks(current_page, total_pages, data_per_page)
			page_package = {
				"current_page":current_page,
				"pageLinks":pageLinks,
				"total_pages":total_pages,
				"data_per_page":data_per_page,
				}
			
			return render(request, 'data_body_template.html', {
				"data_on_page": data_package, #Includes data indexes
				"page_package": page_package, #Includes page links
				"total_data_size": total_data_size,
			})
		
		
		else:
			raise Exception
	except:	
		return HttpResponse("Woopsie!... Something went wrong.")

def data_update(request):
	u = request.user
	if request.method == 'POST' and u.is_authenticated():
		changesMade = simplejson.loads(request.body, "utf-8")
		
		#Get the Lab_Group data to allow direct manipulation.
		saved_data = get_saved_data(u.get_profile().lab_group)	
			
		while (len(changesMade["edit"]) > 0):
			try:
				#An editPackage is [indexChanged, fieldChanged, newValue] 
				editPackage = changesMade["edit"].pop()
				
				indexChanged = int(editPackage[0])-1 #0-based Index
				
				#Translate the DOM class to the Django field.
				fieldChanged = editPackage[1]
				
				#Check that the new value is valid.	
				newValue = editPackage[2] ###CHECK NEW NAMES, CASEY
				assert(quick_validation(fieldChanged, newValue))
				
				#Make the edit in the database
				dataChanged = saved_data[indexChanged]
				setattr(dataChanged, fieldChanged, newValue)
				dataChanged.user = u
				dataChanged.save()
			except:
				print "Invalid data in new edit."
				
		while (len(changesMade["del"]) > 0):###SLOW
			try:
				indexChanged = changesMade["del"].pop()-1 #0-based Index
				saved_data[indexChanged].delete() #Django handles the deletion process.
			except:
				print "ERROR DELETING FILES"
				
		while (len(changesMade["dupl"]) > 0):
			try:
				indexToClone = int(changesMade["dupl"].pop())-1 #0-based Index
				clonedItem = saved_data[indexToClone]
				clonedItem.pk = None #Django creates a new ID for the object.
				clonedItem.user = u #Set the user to the one who performed the duplication.
				clonedItem.save()
			except:
				print "ERROR DUPLICATING FILES"
				
		while (len(changesMade["add"]) > 0): ###Replace "submit" button?
			try:
				indexChanged = changesMade["add"].pop()-1 #0-based Index
				del Lab_Group.saved_data[indexChanged]
			except:
				print "ERROR ADDING FILES"
	return HttpResponse("SUCCESS");

###Register Lab too
def user_registration(request):
	if request.method == "POST":
		uform = UserForm(data = request.POST)
		pform = UserProfileForm(data = request.POST)
		if uform.is_valid() and pform.is_valid():
			#Check that the access_code query for the Lab_Group is correct. 
			q_lab_group = pform.cleaned_data["lab_group"] 
			access_code = Lab_Group.objects.filter(lab_title=q_lab_group)[0].access_code
			
			if pform.cleaned_data["access_code"] == access_code:
				#Create the user to be associated with the profile.
				new_user = uform.save()
				
				#Save the profile
				profile = pform.save(commit = False)
				#Assign the user to the profile
				profile.user = new_user
				profile.save()
				
				new_user = auth.authenticate(username = request.POST["username"],
					password = request.POST["password"])
				auth.login(request, new_user)
				reload_timer = "<div class=\"reload_timer\"></div>"
				return HttpResponse("Registration Successful!"+reload_timer)###Auto-Login?
			else:
				return HttpResponse("Invalid Access Code!")
	else:
		uform = UserForm()
		pform = UserProfileForm()
	return render(request, "registration_form.html", {
		"uform": uform,
		"pform": pform,
	})

def user_login(request):
	if request.method == "POST":
		username = request.POST.get("username", "")
		password = request.POST.get("password", "")
	
		user = auth.authenticate(username=username, password=password)
		
		if user is not None and user.is_active:
			auth.login(request, user)
			#Send the signal to reload the page.
			reload_timer = "<div class=\"reload_timer\"></div>"
			return HttpResponse("Logged in successfully!"+reload_timer)
		else:
			login_fail = True #The login info is not correct.
	else:
		login_fail = False #The user hasn't logged in yet.
	return render(request, "login_form.html", {
		"login_fail": login_fail,
	})
	
def user_logout(request):
	auth.logout(request)
	return HttpResponse("Logged Out!")

#Error Messages:
def display_404_error(request):
	return render(request, '404_error.html')
	
def display_500_error(request):
	return render(request, '500_error.html')
