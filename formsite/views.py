from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.core.cache import cache
from django.shortcuts import render
from django.contrib import auth
from models import *
from validation import *
from djangoappengine.db.utils import get_cursor, set_cursor

import json
import csv
import string
import datetime

import logging###

######################  Controllers  ###################################
class DataManager(object):
	def __init__(self):
		self.data_per_page = 15 #Change to set different numbers of data per page.
		self.current_radius = 4 #Max number of links to display "around" current link.
		self.batch_size = 2
		
	#Collect all data relevant to a specific lab group.
	def collect_all_data(self, lab_group):
		return Data.objects.filter(lab_group=lab_group).order_by("creation_time")

	def calc_total_pages(self, lab_group, total_data_size = None):
		if not total_data_size:
			total_data_size = cache.get("{}|TOTALSIZE".format(lab_group.lab_title))
		total_pages = 1 + int((total_data_size-1)/self.data_per_page)
		if (total_pages < 1): total_pages = 1
		cache.set("{}|TOTALPAGES".format(lab_group.lab_title), total_pages)
		return total_pages
		
	#(Over)write the cursors in the cache.
	def make_cursors(self, lab_group, lab_data = None, page = None, last_page = None):
		if not lab_data.exists():
			lab_data = self.collect_all_data(lab_group)
		if not last_page:
			#Try to retrieve the total number of pages.
			last_page = cache.get("{}|TOTALPAGES".format(lab_group.lab_title))
			if not last_page:
				last_page = self.calc_total_pages(lab_group, lab_data.count())
		if page:
			assert 1 < page <= last_page #Note: The start page should not have a cursor since no elements are skipped.
			assert lab_data.exists() #A page must be given the data starting on that page.
		else:
			page = 2 #Start making cursors from the beginning.
		data_per_cursor = self.data_per_page*self.batch_size
		sub_lab_data = lab_data[0:data_per_cursor]
		logging.info("\n--------Starting while loop!")###
		cursor = None
		while page <= last_page:
			try:
				#If the page doesn't require a cursor, skip to the next page.
				if (page-1) % self.batch_size == 0:
					#Create and cache the cursors.
					cursor = get_cursor(sub_lab_data)
					sub_lab_data = set_cursor(lab_data, cursor)[0:data_per_cursor]
					cache.set("{}|PAGE|CURSOR|{}".format(lab_group.lab_title, page), cursor)
			except Exception as e:
				logging.info("\n---------Oops...\n{}".format(e))###
			page += 1
		if cursor:
			cache.set("{}|PAGE|CURSOR|LAST".format(lab_group.lab_title), cursor)
				
	def find_cursor(self, lab_group, page):
		distance = 0
		cursor = cache.get("{}|PAGE|CURSOR|{}".format(lab_group.lab_title, page))
		while not cursor:
			if page > 1: 
				page -= 1
				distance += 1 #Distance from the page to the closest cursor.
				if distance >= self.batch_size:#If False, cursors need to be made.
					raise Exception("Cursor missing!")
				cursor = cache.get("{}|PAGE|CURSOR|{}".format(lab_group.lab_title, page))
			else:
				cursor="START" #No cursor exists for the page (ie, begin at the "start" of lab_data).
		return cursor, distance	
		
	def get_page_links(self, current, total_pages):
		#Always display the first page.
		page_links = {1} #Use a set to remove any duplicates. 
		
		if total_pages > 1:
			for i in range(current - self.current_radius, current + self.current_radius+1):
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
						i+=1 #Extra addition to account for the new "..." element.
				except:
					pass #At the last element of the list or an error, just skip.
				i += 1
		return page_links		
	
	#Returns the data relevant to a given page (from a cursor).
	#	Note: Will skip lookup if given lab_data or if not set to overwrite.
	def retrieve_data(self, user, page = 1, overwrite = False, lab_data = None):
		if not user.is_authenticated():
			raise Exception("User not logged in.")
		lab_group = user.get_profile().lab_group
		
		logging.info("\n---------Starting retrieve_data!")
		
		#Only retrieve lab_data if it is not already available.
		cached_data = cache.get("{}|PAGEDATA|{}".format(lab_group.lab_title, page))
		
		
		if not cached_data or overwrite:
			if not lab_data.exists():
				logging.info("\n---------Finding cursor!")
				lab_data = self.collect_all_data(user.get_profile().lab_group)
			#Retrieve the closest cursor.
			logging.info("\n---------Finding cursor!")
			cursor, distance = self.find_cursor(lab_group, page)
			#Apply the cursor to the lab data if applicable.###
			logging.info("\n---------Finding BATCH!")
			if cursor!="START":
				batch_lab_data = set_cursor(lab_data, cursor)###
			else:
				batch_lab_data = lab_data
			#Only return the data on the given page.
			logging.info("\n---------Finding Rel Data!")
			rel_lab_data = batch_lab_data[distance*self.data_per_page:(distance+1)*self.data_per_page]
			#Overwrite the existing cache entry.
			cache.set("{}|PAGEDATA|{}".format(lab_group.lab_title, page), list(rel_lab_data))
			logging.info("\n---------Cache set!")
		else: 
			rel_lab_data = cached_data 
			logging.info("\n---------Loaded page from cache!")
		return rel_lab_data
			
	def update_cursors(self, lab_group, index_updated, lab_data=None): #Takes a 0-based index.
		try:
			logging.info("\n---------UPDATING CURSORS!\n")		
			#Find the cursor associated with the updated index:
			page = (index_updated/self.data_per_page) + 1
			if not lab_data.exists():
				lab_data = self.collect_all_data(user.get_profile().lab_group)
			cache.set("{}|TOTALSIZE".format(lab_group.lab_title), lab_data.count())
			if page == 1:
				sub_lab_data = lab_data
			else:
				cursor, distance = self.find_cursor(lab_group, page)
				sub_lab_data = set_cursor(lab_data, cursor)[0:] #Ignore data before the current cursor.
				logging.info("\n---------MEEEEEEEEEEEEEEP...")		
				sub_lab_data = sub_lab_data[distance*self.data_per_page:(distance+1)*self.data_per_page]
			#Make the cursors that follow the current cursor.
				
			#Erase the proceeding cached pages.
			logging.info("\n---------Clearing cash...")		
			total_pages = control.calc_total_pages(lab_group, lab_data.count())
			for i in xrange(page, total_pages+1):
				cache.set("{}|PAGEDATA|{}".format(lab_group.lab_title, i), None)
			logging.info("\n---------Cleared cash on from: {}".format(page))		
			self.make_cursors(lab_group, sub_lab_data, page, total_pages)
		except Exception as e:
			logging.info("\n---------Could not update cursor for page: {}".format(page))		
		
	def get_fresh_page_info(self, request, current_page = None):
		u = request.user
		if not u.is_authenticated():
			raise Exception("User not logged in.")
		try:
			#Gather necessary information from the user's session:
			if not current_page:
				try:
					current_page = int(request.COOKIES.get("current_page"))
				except:
					current_page = 1 #If no page is known, assume 1.
			lab_group = u.get_profile().lab_group
			lab_data = self.collect_all_data(lab_group)
			total_data_size = lab_data.count()
			cache.set("{}|TOTALSIZE".format(lab_group.lab_title), total_data_size)
			total_pages = self.calc_total_pages(lab_group, total_data_size)
			
			#Make sure the page is a valid page.
			if not (0 < current_page <= total_pages):
				current_page = total_pages
			
			#Pack up the session info:
			session = {}
			try:
				logging.info("Trying!")
				session["relevant_data"] = self.retrieve_data(u, current_page, False, lab_data)
				logging.info("Success Loading Old!")
			except:
				#Construct the cursors if necessary.
				logging.info("Making new!")
				self.make_cursors(lab_group, lab_data, None, total_pages)
				session["relevant_data"] = self.retrieve_data(u, current_page, False, lab_data)
				logging.info("Success Loading New!")
			session["total_data_size"] = total_data_size
			session["total_pages"] = total_pages
			session["page_links"] = self.get_page_links(current_page, total_pages)
			session["current_page"] = current_page
			return session
		except Exception as e:
			raise Exception("-Data could not be retrieved for page {}\n--{}.".format(current_page, e))

control = DataManager()
	

######################  Core Views  ####################################
def database(request, control = control):
	#Get all saved data if it exists.
	u = request.user
	try: 
		assert(u.is_authenticated())
		session = control.get_fresh_page_info(request)
		relevant_data = session["relevant_data"]
		total_pages = session["total_pages"]
		page_links = session["page_links"]
		current_page = session["current_page"]
		total_data_size = session["total_data_size"]
	except Exception as e:
		logging.info("\n\nCOULD NOT LOAD SESSION INFO\n{}\n\n".format(e))###
		relevant_data = Data.objects.none() #Don't query the database if U is not logged in.
		total_pages = 1
		page_links = [1]
		current_page = 1
		total_data_size = 0
	
	#Show the overall index of each datum.
	start_index = (current_page-1)*control.data_per_page + 1
	end_index = (current_page)*control.data_per_page + 1
	
	#Prepare packages.
	data_package = zip(relevant_data, range(start_index, end_index))
	page_package = {
		"current_page":current_page,
		"total_pages":total_pages,
		"data_per_page":control.data_per_page,
		"page_links":page_links,
		}
			
	return render(request, 'database_global.html', {
		"data_on_page": data_package, #Includes data and data_indexes.
		"page_package": page_package, 
		"total_data_size": total_data_size,
	})


#Send/receive the compound guide form:
def compound_guide_form(request): #If no data is entered, stay on the current page.
	u = request.user
	success = False
	if u.is_authenticated():
		lab_group = u.get_profile().lab_group
		if request.method == 'POST': 
			#Bind the user's data and verify that it is legit.
			form = CompoundGuideForm(lab_group=lab_group, data=request.POST)
			if form.is_valid():
				#If all data is valid, save the entry.
				form.save()
				#Clear the cache.
				cache.set("{}|COMPOUNDGUIDE".format(lab_group.lab_title), None)
				success = True #Used to display the ribbonMessage.
		else:
			#Submit a blank form if one was not just submitted.
			form = CompoundGuideForm()
	
		guide = collect_CG_entries(u.get_profile().lab_group)
		
		return render(request, 'compound_guide_cell.html', {
			"guide": guide,
			"form": form,
			"success": success,
		})
	else:
		return HttpResponse("Please log in to access the compound guide!")

def edit_CG_entry(request):
	u = request.user
	if request.method == 'POST' and u.is_authenticated():
		changesMade = json.loads(request.body, "utf-8")
		
		#Get the Lab_Group data to allow direct manipulation.
		lab_group = u.get_profile().lab_group
		CG_data = collect_CG_entries(lab_group)	
		
		#Clear the cache.
		cache.set("{}|COMPOUNDGUIDE".format(lab_group.lab_title), None)
		
		#Since only deletions are supported currently. ###
		for index in changesMade:
			try:
				CG_data[int(index)].delete()
				logging.info("\n\nDeleted!")
			except Exception as e:
				logging.info("\n\nCould not delete index! {}".format(e))
			
	return HttpResponse("OK")

#Send/receive the data-entry form:
def data_form(request): #If no data is entered, stay on the current page.
	u = request.user
	success = False
	if request.method == 'POST' and u.is_authenticated(): 
		#Bind the user's data and verify that it is legit.
		form = DataEntryForm(user=u, data=request.POST)
		if form.is_valid():
			#If all data is valid, save the entry.
			form.save()
			cache.set("{}|PAGEDATA|{}".format(lab_group.lab_title, page), None)
			success = True #Used to display the ribbonMessage.
	else:
		#Submit a blank form if one was not just submitted.
		form = DataEntryForm(
			initial={"leak":"No"}
		)
	
	return render(request, 'data_form.html', {
		"form": form,
		"success": success,
	})

#Send/receive the upload-CSV form:
def upload_CSV(request):
	u = request.user
	if request.method == 'POST' and u.is_authenticated(): 
		return upload_data(request)
	else:
		return render(request, 'upload_form.html')
	
######################  Helper Functions ###############################

#Returns a related data entry field (eg, "reactant 1 name" --> "reactant_1")
def get_related_field(heading): ###Not re-read.
	#Strip all punctuation, capitalization, and spacing from the header. 
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
		raise Exception("No valid heading found for \"{}\".".format(heading))###Possible Raise?
	return related_field
		
def live_name_validation(request): ###
	u = request.user
	try:
		if u.is_authenticated() and request.method=="POST":
			saved_data = get_saved_data(u.get_profile().lab_group)
			
			#Give the specific index requested.
			raw_name = json.loads(request.POST["indexRequested"])
			entry = saved_data[index_requested]
			valid_name = ""
		return HttpResponse(valid_name)	
	except:
		return HttpResponse("False")	

######################  Upload/Download Functionality ##################
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
			blacklist = {"x", "-1", -1, "z", "?", "", " "} #Implies absence of data.
			unknown_label = "?" #The label that blacklist values will inherit.
			not_required_label = "" #The label that auto-added values will inherit.
			
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
						recheck_fields = set()
						for i in range(1,6): #Since only 5 reactants supported...
							if not "unit_{}".format(i) in set_user_fields:
								user_fields.insert(
									user_fields.index("quantity_{}".format(i))+1,
										"unit_{}".format(i))
								auto_added_fields.add("unit_{}".format(i))
								
						#Assert that there are no duplicates in the list
						#	and that all fields are valid.
						if len(user_fields) != len(true_fields):
							raise Exception("Too many columns present! {} needed but {} found.".format(len(true_fields), len(user_fields)))
							
						for field in true_fields:
							assert(field in user_fields)
							
						headings_valid = True
						continue
				except Exception as e:
					validation_attempt += 1
					error_log += "<li>{}</li>".format(e)
					continue
					
				#All other rows should have data corresponding to the headings.
				try:
					#Create new object that will receive the uploaded data.
					entry_fields = {}
					#Access the data from the uploaded file.
					i = 0 # data_group index (gets a datum)
					j = 0 # user_fields index. (gets a field)
					#
					while (i < len(data_group) or j < len(user_fields)):
						#Required since data and fields may be disjunct from missing units.
						datum = data_group[i]
						field = user_fields[j]
						if field in auto_added_fields: 
							#Skip the field if it was auto-added because
							#	no data is present in the generated column.
							j += 1
							try:
								entry_fields[field] #Check if the field exists already.
							except:
								entry_fields[field] = not_required_label
							continue
						try:
							#If the datum isn't helpful, don't remember it.
							if datum in blacklist: 
								if field in not_required:
									datum = not_required_label
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
					create_data_entry(u, **entry_fields)
					added_quantity += 1
				except Exception as e:
					error_log += "<li>{} ---- Failed assigning \"{}\" as \"{}\".</li>".format(e, datum, field)
					error_quantity +=1
				
			#Update cursors if data was added.
			if added_quantity:
				lab_group = u.get_profile().lab_group
				lab_data = control.collect_all_data(lab_group)
				first_index_updated = lab_data.count()-added_quantity
				control.update_cursors(lab_group, first_index_updated, lab_data)
		except Exception as e:
			error_log += "<li>{}</li>".format(e)
	
	###Bulk Upload Instead?
	###Template?
		message = "<h1>Results:</h1>".format(added_quantity)
		message += "Added: {}<br/>".format(added_quantity)
		message += "Failed: {}".format(error_quantity)
		message += "<h1>Error Log: </h1>{}".format(error_log)
		message += "<a href=\"/database/\">Return to Database</a>"
		return HttpResponse(message)
	elif not u.is_authenticated():
		return HttpResponse("<p>Please log in to upload data.</p>")
	else:
		return HttpResponse("<p>No file uploaded.</p>")

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
		lab_data = control.collect_all_data(u.get_profile().lab_group)
		
		for entry in lab_data:
			row = []
			try:
				for field in headers:
					row += [eval("entry.{}".format(field))]
				writer.writerow(row)
			except:
				pass
		return CSV_file #ie, return HttpResponse(content_type="text/csv")
	else:
		return HttpResponse("<p>Please log in to download data.</p>")

######################  Change Page ####################################
def data_transmit(request, num = 1, control=control):
	try:
		if request.method == "GET" and request.user.is_authenticated():
			#Get the request information.
			u = request.user
			session = control.get_fresh_page_info(request, int(num))
			
			#Get the necessary data for a page change.
			relevant_data = session["relevant_data"]
			total_pages = session["total_pages"]
			page_links = session["page_links"]
			current_page = session["current_page"]
			total_data_size = session["total_data_size"]
			
			#Only send the data on the requested page. Note that Lab_Group.saved_data is a 0-based index.
			start_index = (current_page-1)*control.data_per_page 
			end_index = (current_page)*control.data_per_page
			index_range = range(start_index+1, end_index+1) #1-based index for IDs (para los users).
		
			#Prepare packages.
			data_package = zip(relevant_data, index_range)
			page_package = {
				"current_page":current_page,
				"total_pages":total_pages,
				"data_per_page":control.data_per_page,
				"page_links":page_links,
				}
				
			return render(request, 'data_and_page_container.html', {
				"data_on_page": data_package, #Includes data indexes
				"page_package": page_package, #Includes page links
				"total_data_size": total_data_size,
			})
		else:
			return HttpResponse("<p>Please log in to view your data.</p>")
	except:	
		return HttpResponse("<p>Woopsie!... Something went wrong.</p>")

######################  Data Transmit ##################################
#Send the CG name pairs to the client.
def send_CG_names(request):
	u = request.user
	if u.is_authenticated():
		name_pairs = collect_CG_name_pairs(lab_group, overwrite=False)
		return HttpResponse(name_pairs)
	return HttpResponse("Please log in to see data.")

######################  Update Data ####################################
def data_update(request, control=control): ###Lump together?
	u = request.user
	if request.method == 'POST' and u.is_authenticated():
		changesMade = json.loads(request.body, "utf-8")
		
		#Get the Lab_Group data to allow direct manipulation.
		lab_group = u.get_profile().lab_group
		lab_data = control.collect_all_data(lab_group)
		total_size = lab_data.count()	
		earliest_index = -1 #Impossible, so if it changes, then a change exists.
			
		while (len(changesMade["edit"]) > 0):
			try:
				#An editPackage is [indexChanged, fieldChanged, newValue]. 
				editPackage = changesMade["edit"].pop()
				indexChanged = int(editPackage[0])-1 #Translate to 0-based Index
				fieldChanged = editPackage[1] 
				newValue = editPackage[2] #Edits are validated client-side.
				
				#Make the edit in the database
				dataChanged = lab_data[indexChanged]
				setattr(dataChanged, fieldChanged, newValue)
				dataChanged.user = u
				dataChanged.save()
				
				#Erase the cached page since it is no longer up-to-date.
				page = (indexChanged/control.data_per_page) + 1
				cache.set("{}|PAGEDATA|{}".format(lab_group.lab_title, page), None)
			except:
				pass
		while (len(changesMade["dupl"]) > 0):
			try:
				indexToClone = int(changesMade["dupl"].pop())-1 #0-based Index
				clonedItem = lab_data[indexToClone]
				clonedItem.pk = None #Django creates a new ID for the object.
				clonedItem.user = u #Set the user to the one who performed the duplication.
				clonedItem.creation_time = str(datetime.datetime.now())
				clonedItem.save()
				if earliest_index < total_size or earliest_index == -1: 
					earliest_index = total_size
				total_size += 1
			except:
				pass
		while (len(changesMade["del"]) > 0):###SLOW
			try:
				indexChanged = changesMade["del"].pop()-1 #0-based Index
				lab_data[indexChanged].delete() #Django handles the deletion process.
				if earliest_index < indexChanged or earliest_index == -1: 
					earliest_index = indexChanged
			except:
				pass
		if earliest_index != -1:
			control.update_cursors(lab_group, earliest_index, lab_data)
		
	return HttpResponse("OK"); #Django requires an HttpResponse...

def get_full_datum(request, control=control):
	u = request.user
	try:
		if u.is_authenticated() and request.method=="POST":
			#Give the specific index requested.
			index_requested = json.loads(request.POST["indexRequested"])
			###Is this efficient?
			entry = control.collect_all_data(u.get_profile().lab_group)[index_requested]
			
		return render(request, "full_datum.html", {
			"entry":entry
		})
	except Exception as e:
		return HttpResponse("<p>Data could not be loaded!</p>")		

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
	#"render" is needed to properly display template:
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
	
def lab_registration(request): ###Not finished.
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
