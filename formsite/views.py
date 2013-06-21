from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.shortcuts import render
from django.contrib import auth
from models import *
from validation import *

import json
import csv
import string
import datetime

######################  Static  ########################################
class Controller(object):
	def __init__(self):
		self.data_per_page = 15 #Change to set different numbers of data per page.
		self.current_radius = 4 #Max number of links to display "around" current link.
	
control = Controller()	

######################  Session Info  ####################################
#Returns all data that belongs to a specific lab group.
def get_saved_data(lab_group):
	return Data.objects.filter(lab_group=lab_group)

def get_total_pages(total_data_size):
	total_pages = int(1 + (total_data_size-1)/control.data_per_page)
	if (total_pages < 1): total_pages = 1
	return total_pages
		
def get_page_links(current, total_pages):
	#Always display the first page.
	page_links = {1} #Use a set to remove any duplicates. 
	
	if total_pages > 1:
		for i in range(current - control.current_radius, current + control.current_radius+1):
			if (1 < i < total_pages): page_links.add(i)
			
		#Always display the last page if applicable.
		page_links.add(total_pages)
		
	#Convert page_links to an ordered list.
	page_links = list(page_links)
	page_links.sort()
	
	if len(page_links) >= 2:
		i=0
		while i < len(page_links):
			try:
				if page_links[i+1]-page_links[i] > 1: #If a gap exists between two numbers, add an ellipsis. 
					page_links.insert(i+1,"...")
					i+=1 #Extra addition to account for new "..." element.
			except:
				pass #At the last element of the list or an error, just skip.
			i += 1
	return page_links
	
#Used to get dump of up-to-date session information. 
def get_fresh_session_info(current, lab_group):
	session = {}
	session["saved_data"] = get_saved_data(lab_group)
	session["total_pages"] = get_total_pages(len(session["saved_data"]))
	session["page_links"] = get_page_links(current, session["total_pages"])
	return session
	

######################  Core Views  ####################################
def database(request, control = control):
	#Get all saved data if it exists.
	if request.user.is_authenticated():
		u = request.user
		
		#Get the user's current page.
		try:
			current_page = int(request.COOKIES.get("current_page"))
		except:
			current_page = 1 #Start on the last/latest page.
		
		session = get_fresh_session_info(current_page, u.get_profile().lab_group)
		saved_data = session["saved_data"]
		total_pages = session["total_pages"]
		page_links = session["page_links"]
		
	else:
		saved_data = [] #Don't query the database if U is not logged in.
		total_pages = 1
		page_links = [1]
		current_page = 1
	
	#Only send the data on the requested page. Note that Lab_Group.saved_data is a 0-based index.
	start_index = (current_page-1)*control.data_per_page 
	end_index = (current_page)*control.data_per_page
	index_range = range(start_index+1, end_index+1) #1-based index for visuals

	#Prepare packages.
	data_package = zip(saved_data[start_index:end_index], index_range)
	page_package = {
		"current_page":current_page,
		"total_pages":total_pages,
		"data_per_page":control.data_per_page,
		"page_links":page_links,
		}
			
	return render(request, 'database_global.html', {
		"data_on_page": data_package, #Includes data and data_indexes.
		"page_package": page_package, 
		"total_data_size": len(saved_data),
	})


#Send/receive the data-entry form:
def data_form(request): #If no data is entered, stay on the current page.
	u = request.user
	success = False
	if request.method == 'POST' and u.is_authenticated(): 
		#Bind the user's data and verify that it is legit.
		form = DataEntryForm(user=u, data=request.POST)
		if form.is_valid():
			#If all data is valid, save the entry (and the submitter)
			form.save()
			success = True
	else:
		#Submit a blank form if one was not just submitted.
		form = DataEntryForm()
	
	return render(request, 'data_form.html', {
		"form": form,
		"success": success,
	})

#Send/receive the upload-CSV form:
def upload_CSV(request):
	u = request.user
	print "FILE:",request.FILES###
	if request.method == 'POST' and u.is_authenticated(): 
		return upload_data(request)
	else:
		return render(request, 'upload_form.html')
	
######################  Helper Functions ###############################

#Helper function that returns a related data entry field ("reactant 1 name" --> "reactant_1")
def get_related_field(heading): ###Not re-read.
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
		
def upload_data(request): ###Not re-read.
	true_fields = get_data_field_names()
	not_required = { ###Auto-generate?
			"reactant_3", "quantity_3", "unit_3", 
			"reactant_4", "quantity_4", "unit_4",
			"reactant_5", "quantity_5", "unit_5",
			"notes"
		}
			
	u = request.user
	
	if request.method=="POST" and request.FILES and u.is_authenticated():
		uploaded_file = request.FILES["file"]
		added_quantity=0
		error_quantity=0
		error_log = ""
		
		try:
			blacklist = {"x", "-1", -1, "z", "?", "", " "} #Implies absence of data. ###
			unknown_label = "?" #The label that blacklist values will inherit.
			
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
							if field[:-2]=="quantity":
								auto_add
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
									if field == "quantity_5": print "YES!"
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
					print "1"###
					print entry_fields
					print "<br />"
					create_data_entry(u, **entry_fields)
					added_quantity += 1###
				except Exception as e:
					error_log += "ERROR:",e,"\n",data_group,"at",field,"({})<br/>".format(datum)
					error_quantity +=1###
					
		except Exception as e:
			print "Data could not be uploaded."
	
	###Bulk Upload Instead?
		message = "<h1>Results:</h1>".format(added_quantity)
		message += "Added: {}<br/>".format(added_quantity)
		message += "Failed: {}".format(error_quantity)
		message += "<h1>Error Log: </h1>{}".format(error_log)
		return HttpResponse(message)
	else:
		return HttpResponse("<p>Please log in to upload data.</p>")

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
	else:
		return HttpResponse("<p>Please log in to download data.</p>")
	

######################  Change Page ####################################
def data_transmit(request, num = 0, control=control):
	try:
		if request.method == "GET" and request.user.is_authenticated():
			#Get the request information.
			u = request.user
			requested_page = int(num)
				
			#Get the necessary data for a page change.
			session = get_fresh_session_info(requested_page, u.get_profile().lab_group)
			saved_data = session["saved_data"]
			total_pages = session["total_pages"]
			page_links = session["page_links"]
			
			#If the page does not exist, raise a 404. 
			try:
				assert 0 < requested_page <= total_pages
			except:
				raise Http404
		
			#Only send the data on the requested page. Note that Lab_Group.saved_data is a 0-based index.
			start_index = (requested_page-1)*control.data_per_page 
			end_index = (requested_page)*control.data_per_page
			index_range = range(start_index+1, end_index+1) #1-based index for visuals
		
			#Prepare packages.
			data_package = zip(saved_data[start_index:end_index], index_range)
			page_package = {
				"current_page":requested_page,
				"total_pages":total_pages,
				"data_per_page":control.data_per_page,
				"page_links":page_links,
				}
				
			return render(request, 'data_and_page_container.html', {
				"data_on_page": data_package, #Includes data indexes
				"page_package": page_package, #Includes page links
				"total_data_size": len(saved_data),
			})
		else:
			return HttpResponse("<p>Please log in to view your data.</p>")
	except:	
		return HttpResponse("<p>Woopsie!... Something went wrong.</p>")

######################  Update Data ####################################
def data_update(request):
	u = request.user
	if request.method == 'POST' and u.is_authenticated():
		changesMade = json.loads(request.body, "utf-8")
		
		#Get the Lab_Group data to allow direct manipulation.
		saved_data = get_saved_data(u.get_profile().lab_group)	
			
		while (len(changesMade["edit"]) > 0):
			try:
				#An editPackage is [indexChanged, fieldChanged, newValue]. 
				editPackage = changesMade["edit"].pop()
				indexChanged = int(editPackage[0])-1 #Translate to 0-based Index
				fieldChanged = editPackage[1] 
				newValue = editPackage[2] ###CHECK NEW NAMES, CASEY
				assert(quick_validation(fieldChanged, newValue)) #Check that the new value is valid.
				
				#Make the edit in the database
				dataChanged = saved_data[indexChanged]
				setattr(dataChanged, fieldChanged, newValue)
				dataChanged.user = u
				dataChanged.save()
			except:
				pass
		while (len(changesMade["del"]) > 0):###SLOW
			try:
				indexChanged = changesMade["del"].pop()-1 #0-based Index
				saved_data[indexChanged].delete() #Django handles the deletion process.
			except:
				pass
		while (len(changesMade["dupl"]) > 0):
			try:
				indexToClone = int(changesMade["dupl"].pop())-1 #0-based Index
				clonedItem = saved_data[indexToClone]
				clonedItem.pk = None #Django creates a new ID for the object.
				clonedItem.user = u #Set the user to the one who performed the duplication.
				clonedItem.save()
			except:
				pass
	return HttpResponse("SUCCESS"); #Django requires an HttpResponse...

######################  User Auth ######################################
def user_login(request):
	login_fail = False #The user hasn't logged in yet...
	
	if request.method == "POST":
		username = request.POST.get("username", "")
		password = request.POST.get("password", "")
		user = auth.authenticate(username=username, password=password)
		if user is not None and user.is_active:
			auth.login(request, user)
			
			return HttpResponse("Logged in successfully! <div class=reloadActivator></div>"); #Only the reloadActivator is "required" here.
		else:
			login_fail = True #The login info is not correct.
	return render(request, "login_form.html", {
		"login_fail": login_fail,
	})
	
def user_logout(request):
	auth.logout(request)
	return HttpResponse("Logged Out!")

#Redirects user to the appropriate registration screen.
def registration_prompt(request):
	#Needed to properly render template.
	return render(request, "registration_cell.html", {})

def user_registration(request):
	if request.method == "POST":
		form = [UserForm(data = request.POST), UserProfileForm(data = request.POST)]
		if form[0].is_valid() and form[1].is_valid():
			#Check that the access_code query for the Lab_Group is correct. 
			lab_group = form[1].cleaned_data["lab_group"] 
			access_code = Lab_Group.objects.filter(lab_title=lab_group)[0].access_code
			
			if form[1].cleaned_data["access_code"] == access_code:
				#Create the user to be associated with the profile.
				new_user = form[0].save()
				#Save the profile
				profile = form[1].save(commit = False)
				#Assign the user to the profile
				profile.user = new_user
				profile.save()
				#Politely log the user in!
				new_user = auth.authenticate(username = request.POST["username"],
					password = request.POST["password"])
				auth.login(request, new_user)
				reload_timer = "<div class=reloadActivator></div>"
				return HttpResponse("Registration Successful!"+reload_timer)
			else:
				return HttpResponse("Invalid Access Code!")
	else:
		form = [UserForm(), UserProfileForm()]
	return render(request, "user_registration_form.html", {
		"user_form": form[0],
		"profile_form": form[1],
	})
	
def lab_registration(request):
	if request.method=="POST":
		pass
	else:
		####form = LabForm()
		pass
	return render(request, "lab_registration_form.html", {
		###"lab_form": form,
	})	

######################  Error Messages  ################################
def display_404_error(request):
	return render(request, '404_error.html')
	
def display_500_error(request):
	return render(request, '500_error.html')
