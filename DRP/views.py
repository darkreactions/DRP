from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.shortcuts import render
from django.contrib import auth
from django.db.models import Q
from models import *
from validation import *

import json
import csv
import string
import datetime

from svg_construction import *
#from construct_descriptor_table import *
from data_config import CONFIG


######################  Controllers  ###################################

#Manages pages and data:
class DataManager(object):
	def __init__(self):
		self.data_per_page = 15 #Change to set different numbers of data per page.
		self.current_radius = 4 #Max number of links to display "around" current link.

	#Collect all data relevant to a specific lab group.
	def collect_all_data(self, lab_group):
		return Data.objects.filter(lab_group=lab_group).order_by("creation_time")

	def calc_total_pages(self, lab_group, total_data_size = None):
		if not total_data_size:
			total_data_size = get_cache(lab_group, "TOTALSIZE")
		total_pages = 1 + int((total_data_size-1)/self.data_per_page)
		if (total_pages < 1): total_pages = 1
		set_cache(lab_group, "TOTALPAGES", total_pages)
		return total_pages

	#Erases the cached page which contains a specific index.
	def clear_page_of(self, lab_group, indexChanged):
		page = (indexChanged/self.data_per_page) + 1
		set_cache(lab_group, "PAGEDATA|{}".format(page), None)

	def clear_all_page_caches(self, lab_group):
		for page in xrange(1, self.calc_total_pages(lab_group)+1):
			set_cache(lab_group, "PAGEDATA|{}".format(page), None)

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

		#Only retrieve lab_data if it is not already available.
		cached_data = get_cache(lab_group, "PAGEDATA|{}".format(page))

		if not cached_data or overwrite:
			if not lab_data.exists():
				lab_data = self.collect_all_data(user.get_profile().lab_group)
			#Apply the cursor to the lab data if applicable.###
			#Only return the data on the given page.
			rel_lab_data = lab_data[(page-1)*self.data_per_page:(page)*self.data_per_page]
			#Overwrite the existing cache entry.
			set_cache(lab_group, "PAGEDATA|{}".format(page), list(rel_lab_data))
		else:
			rel_lab_data = cached_data
		return rel_lab_data

	def get_page_info(self, request, current_page = None):
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
			set_cache(lab_group, "TOTALSIZE", total_data_size)
			total_pages = self.calc_total_pages(lab_group, total_data_size)

			#Make sure the page is a valid page.
			if not (0 < current_page <= total_pages):
				current_page = total_pages

			#Pack up the session info:
			session = {
				"relevant_data": self.retrieve_data(u, current_page, False, lab_data),
				"total_data_size": total_data_size,
				"total_pages": total_pages,
				"page_links": self.get_page_links(current_page, total_pages),
				"current_page": current_page,
			}
			return session
		except Exception as e:
			raise Exception("-Data could not be retrieved for page {}\n--{}.".format(current_page, e))

def revalidate_all_data(lab_group, invalid_only = True):
	data_to_validate = Data.objects.filter(lab_group=lab_group)

	if invalid_only:
		data_to_validate = data_to_validate.filter(is_valid=False)

	print "Found {} to validate.".format(data_to_validate.count())

	if data_to_validate.exists():
		for data in data_to_validate:
			revalidate_data(data, lab_group, batch=True)
		#Clear the page caches
		control.clear_all_page_caches(lab_group)
control = DataManager()

######################  Analysis Functions  ####################################
#Return all data before/after a specific date (ignoring time).
def filter_by_date(lab_data, raw_date, direction="after"):
	#Convert the date input into a usable string.
	#Note, date must be given as MM-DD-YY.
	date = datetime.datetime.strptime(raw_date, "%m-%d-%y")
	date_string = str(date)

	#Get the reactions before/after a specific date.
	if direction.lower() == "after":
		print "Getting data AFTER {}".format(date_string)###
		filtered_data = lab_data.filter(creation_time__gte=date_string)
	else:
		print "Getting data BEFORE {}".format(date_string)###
		filtered_data = lab_data.filter(creation_time__lte=date_string)

	return filtered_data


######################  Core Views  ####################################
def database(request, control = control):
	#Get user data if it exists.
	u = request.user
	if u.is_authenticated():
		session = control.get_page_info(request)
		relevant_data = session["relevant_data"]
		total_pages = session["total_pages"]
		page_links = session["page_links"]
		current_page = session["current_page"]
		total_data_size = session["total_data_size"]
	else:
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

###
import time

def predictions(request):
	#Variable Setup
	u = request.user
	fatal_message = ""
	svg = ""

	if u.is_authenticated():
		try:
			lab_group = u.get_profile().lab_group
			svg = get_cache(lab_group, "TESTSVG")###
			if not svg:
				start_time = time.time()###
				#Attempt to validate any invalid data.
				###revalidate_all_data(lab_group) ###Validates all data or just user data?

				#Create and cache the SVG.
				svg = generate_svg(u.get_profile().lab_group)
				print "Took {} seconds overall.".format(time.time()-start_time)###

				set_cache(lab_group, "TESTSVG", svg, 86400)
		except Exception as e:
			fatal_message = e
		#construct_descriptor_table("cat","dog")
	else:
		fatal_message = "Please log in to view predictions."

	return render(request, 'predictions_global.html', {
		"fatal_message": fatal_message,
		"svg": svg, #Includes data and data_indexes.
	})

def gather_SVG(request):
	u = request.user
	fatal_message = ""
	svg = ""

	if u.is_authenticated() and request.method=="POST":
		try:
			step = request.POST.get("step")
			source = request.POST.get("source")
			svg = generate_svg(u.get_profile().lab_group, step, source)
		except Exception as e:
			fatal_message = e
		#construct_descriptor_table("cat","dog")
	elif request.method!="POST":
		fatal_message = "Could not gather information about SVG to create."
	else:
		fatal_message = "Please log in to view predictions."

	if fatal_message:
		return HttpResponse("<p class=\"fatalError\">{}</p>".format(fatal_message))
	return HttpResponse(svg)


######################  Searching  #####################################
def search(request):
	u = request.user
	max_search_size = 1000
	search_page_size = 150
	if u.is_authenticated():
		other_fields = get_model_field_names(both=True, unique_only=True)
		if request.method=="POST":

			query_list = json.loads(request.POST.get("current_query"))
			if not query_list:
				return HttpResponse("No search specified!")

			lab_group = u.get_profile().lab_group
			all_data = Data.objects.filter(lab_group=lab_group) #Order doesn't matter at this point in search.
			filters = ""

			for query in query_list:
				field = query.get("field")
				value = query.get("value")

				#Check the field and value for possible attacks. ###Make better before you show Paul or he will cry.
				for char in field:
					assert(not char in "\"'!={}/\\~`")
				for char in value:
					assert(not char in "\"'!={}/\\~`")

				if field in list_fields:
					Q_string = ""
					for i in CONFIG.reactant_range():
						Q_string += "Q({}_{}=\"{}\")|".format(field, i, value)
					Q_string = Q_string[:-1] #Remove the last trailing "|".
					filters += ".filter({})".format(Q_string)
				else:
					filters += ".filter({}=\"{}\")".format(field, value)

			entries = eval("all_data{}".format(filters)).order_by("creation_time")[:max_search_size]
			pk_list = entries[:max_search_size].values("pk")

			return render(request, 'search_results.html', {
				"entries": entries[:search_page_size], #Only give back a specific amount of entries.
				"search_page_size":search_page_size,
				"max_search_size": max_search_size,
				"entries_found": entries.count(),
				"no_search": False,
				"other_fields": other_fields,
			})

		else:
			return render(request, 'search_global.html', {
			"no_search": True,
			"other_fields": other_fields,
			})
	else:
		return HttpResponse("Log in to search data.")

######################  CG Guide  ######################################
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
				#Clear the cached CG data.
				set_cache(lab_group, "COMPOUNDGUIDE", None)
				set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", None)
				success = True #Used to display the ribbonMessage.
		else:
			#Submit a blank form if one was not just submitted.
			form = CompoundGuideForm()

		guide = collect_CG_entries(lab_group)

		return render(request, 'compound_guide_cell.html', {
			"guide": guide,
			"form": form,
			"success": success,
		})
	else:
		return HttpResponse("Please log in to access the compound guide!")

def edit_CG_entry(request): ###Edits?
	u = request.user
	if request.method == 'POST' and u.is_authenticated():
		changesMade = json.loads(request.body, "utf-8")

		#Get the Lab_Group data to allow direct manipulation.
		lab_group = u.get_profile().lab_group
		CG_data = collect_CG_entries(lab_group)

		#Clear the cached CG entries.
		set_cache(lab_group, "COMPOUNDGUIDE", None)
		set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", None)

		###Should make cleaner and "directed" rather than "catching", Casey --Tired Casey
		edit_type = changesMade.get("type", "del")

		if changesMade["type"]=="del":
			#Since only deletions are supported currently. ###
			for index in changesMade: ###Tie names to CG abbrevs directly?
				try:
					CG_data[int(index)].delete()
				except Exception as e:
					print("\n\nCould not delete index: {}".format(e))
		elif changesMade["type"]=="edit":
			try:
				#Variable Setup
				field = changesMade["field"]
				new_val  = changesMade["newVal"]
				old_val  = changesMade["oldVal"]
				compound = changesMade["compound"]

				assert(new_val and compound and field)

				#Gather all relevant data (and check for duplicates).
				try:
					changed_entry = CompoundEntry.objects.get(lab_group=lab_group, compound=compound)
				except:
					return HttpResponse("Please delete duplicate entry!")

				#Make sure the compound isn't already being used.
				if field=="compound":
					possible_entry = CompoundEntry.objects.filter(lab_group=lab_group, compound=new_val)
					if possible_entry.exists():
						return HttpResponse("Already used!")
				elif field=="abbrev":
					possible_entry = CompoundEntry.objects.filter(lab_group=lab_group, abbrev=new_val)
					if possible_entry.exists():
						return HttpResponse("Already used!")


				#Make sure the new value does not invalidate the entry.
				dirty_data = model_to_dict(changed_entry)
				dirty_data[field] = new_val
				clean_data, errors = CG_validation(dirty_data, lab_group, editing_this=True)
				new_val = clean_data[field]
				if errors:
					return HttpResponse(errors.values())

				#Commit the change to the Entry.
				setattr(changed_entry, field, new_val)
				changed_entry.save()

				#Make any edits to the Data if needed.
				if field=="abbrev":
					#Don't overwrite ALL empty entries -- only those related to the compound.
					if not old_val:
						old_val = compound

					for i in CONFIG.reactant_range():
						exec("old_data = Data.objects.filter(lab_group=lab_group, reactant_{}=old_val)".format(i))
						exec("old_data.update(reactant_{}=new_val)".format(i))

				#Clear the cached CG
				set_cache(lab_group, "COMPOUNDGUIDE", None)
				set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", None)
			except Exception as e:
				print e###
				return HttpResponse("Invalid!")

	return HttpResponse("0")

	##################  Helper Functions ###############################
def guess_type(datum):
	guess = ""
	datum=datum.lower()
	if "wat" in datum or "h2o" in datum:
		return "Water"
	if "oxa" in datum:
		return "Ox"
	if ("eth" in datum or "prop" in datum or "but" in datum or "amin" in datum
		or "pip" in datum or ("c" in datum and not "cl" in datum)):
		return "Org"
	if "ol" in datum:
		return "Sol"
	return "Inorg" #Default to inorganic if no guess is uncovered.

######################  Database Functions  ############################
#Send/receive the data-entry form:
def data_form(request, copy_index=None): #If no data is entered, stay on the current page.
	u = request.user
	success = False
	if request.method == 'POST' and u.is_authenticated():
		#Bind the user's data and verify that it is legit.
		form = DataEntryForm(user=u, data=request.POST)
		if form.is_valid():
			#If all data is valid, save the entry.
			form.save()
			lab_group = u.get_profile().lab_group

			#Clear the cache of the last page.
			old_data_size = get_cache(lab_group, "TOTALSIZE")
			control.clear_page_of(lab_group, old_data_size)

			#Refresh the TOTALSIZE cache.
			set_cache(lab_group, "TOTALSIZE", old_data_size + 1)
			set_cache(lab_group, "REFS", None)
			success = True #Used to display the ribbonMessage.
	else:
		#Submit a blank/auto-filled form if one was not just submitted.
		try:
			#Verify that an index is given is present.
			assert(copy_index != None)
			datum = control.collect_all_data(u.get_profile().lab_group)[int(copy_index)-1]
			#Get the data from a specific index.
			initial_fields = {}
			for field in get_model_field_names(model="Data"):
				initial_fields[field] = getattr(datum, field)
			form = DataEntryForm(
				initial=initial_fields
			)

		except:
			form = DataEntryForm(
				initial={"leak":"No"}
			)

	return render(request, 'data_form.html', {
		"form": form,
		"success": success,
	})

	##################  Helper Functions ###############################

#Returns a related data entry field (eg, "reactant 1 name" --> "reactant_1")
def get_related_field(heading, model="Data"): ###Not re-read.
	#Strip all punctuation, capitalization, and spacing from the header.
	stripped_heading = heading.translate(None, string.punctuation)
	stripped_heading = stripped_heading.translate(None, " ").lower()
	stripped_heading = stripped_heading[:20] #Limit the checked heading (saves time if super long).

	if model=="Data":
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
			raise Exception("Not a valid heading: <div class=failedUploadData>{}</div>".format(heading))
	elif model=="CompoundEntry":
		if "cas" in stripped_heading:
			related_field = "CAS_ID"
		elif "type" in stripped_heading:
			related_field = "compound_type"
		elif ("comp" in stripped_heading or "full" in stripped_heading
			or "name" in stripped_heading):
			related_field = "compound"
		elif "abbr" in stripped_heading or "short" in stripped_heading:
			related_field = "abbrev"
		else:
			raise Exception("Not a valid heading: <div class=failedUploadData>{}</div>".format(heading))
	else:
		raise Exception("Unknown model specification for relations.")
	return related_field

######################  Upload/Download   ##############################
def upload_CSV(request, model="Data"): ###Not re-read.
	u = request.user

	#Variable setup.
	error_log = []
	fatal_message = ""
	added_quantity = 0
	error_quantity = 0
	entry_list = [] #Used to bulk create entries.
	no_abbrev = False

	if request.method=="POST" and request.FILES and u.is_authenticated():
		#Get the file and model specification from the POST request.
		lab_group = u.get_profile().lab_group
		model=request.POST["dataType"]
		uploaded_file = request.FILES["file"]

		#Get the appropriate field names.
		if model=="Data":
			not_required = {
					"reactant_3", "quantity_3", "unit_3",
					"reactant_4", "quantity_4", "unit_4",
					"reactant_5", "quantity_5", "unit_5",
					"notes", "duplicate_of", "recommended",
				}
			allow_unknowns = True
		elif model=="CompoundEntry":
			not_required = {
					"CAS_ID"
				}
			allow_unknowns = False
			abbrev_dict = collect_CG_name_pairs(lab_group)
		else:
			raise Exception("Unknown model specified in upload.")

		true_fields = get_model_field_names(model=model)
		row_num = 2 #Assuming row 1 is headings.

		try:
			#Attempt to validate the headings of the uploaded CSV.
			headings_valid = False
			validation_attempt = 0
			#Separate data into groups of fields -- then separate fields.
			for row in csv.reader(uploaded_file, delimiter=","):
				#Variable Setup for each data group.
				data_is_valid = True

				#The first row should be a series of headings.
				if validation_attempt > 5:
					fatal_message = "File doesn't have valid headings."
					raise Exception("Exceeded validation attempts.")
				try:
					if not headings_valid: #Remember which column has which heading.
						#Translate the user's headings into usable field names.
						user_fields = []
						row_length = 0 #Remember how many elements to read per row.
						for field in row:
							if field == "":
								break #Stop checking for headings if a "" is encountered.
							user_fields.append(get_related_field(field, model=model))
							row_length+=1

						#If unit columns were not supplied, add them after each quantity.
						set_user_fields = set(user_fields)
						auto_added_fields = set()
						if model=="Data":
							for i in CONFIG.reactant_range(): #Since only 5 reactants supported...
								if not "unit_{}".format(i) in set_user_fields:
									user_fields.insert(
										user_fields.index("quantity_{}".format(i))+1,
											"unit_{}".format(i))
									auto_added_fields.add("unit_{}".format(i))
							for opt in ["duplicate_of", "recommended"]:
								if not opt in set_user_fields:
									user_fields.append(opt)
									auto_added_fields.add(opt)

						elif model=="CompoundEntry":
							for opt in ["CAS_ID", "abbrev"]:
								if not opt in set_user_fields:
									user_fields.append(opt)
									auto_added_fields.add(opt)

						#Assert that there are no duplicates in the list
						#	and that all fields are valid.
						if len(user_fields) != len(true_fields):
							raise Exception("Invalid column quantity! {} needed but {} found.".format(len(true_fields), len(user_fields)))

						for field in true_fields:
							assert(field in user_fields)
						headings_valid = True
						continue
				except Exception as e:
					validation_attempt += 1
					error_log.append([str(e)])
					continue

				#All other rows should have data corresponding to the headings.
				try:
					#Create new object that will receive the uploaded data.
					model_fields = {}
					#Access the data from the uploaded file.
					i = 0 # row index (gets a datum)
					j = 0 # user_fields index. (gets a field)
					no_abbrev = False

					#Get the set of references so no refs are reused.
					ref_set = get_ref_set(lab_group)

					###while (i < row_length or j < len(user_fields)):
					while (j < len(user_fields)):
						#Required since data and fields may be disjunct from missing units.
						field = user_fields[j]

						#Only get the datum if it is available.
						if field in auto_added_fields:
							#Skip the field if it was auto-added because
							#	no data is present for the generated column.
							j += 1
							try:
								if model=="Data":
									if field=="recommended":
										#Assume blank abbrevs should be same as compound.
										model_fields[field] = "No"
									elif field=="duplicate_of":
										#Assume blank abbrevs should be same as compound.
										model_fields[field] = ""
								if model=="CompoundEntry":
									if field=="abbrev":
										#Assume blank abbrevs should be same as compound.
										model_fields[field] = model_fields["compound"]
								model_fields[field] #Check if the field exists already.
							except:
								model_fields[field] = CONFIG.not_required_label
							continue
						else:
							datum = row[i]

						try:
							#If the datum isn't helpful, don't remember it.
							if type(datum)==str:
								datum = datum.replace("\n","").replace("\t","")#TODO: Find better way?
							if datum in CONFIG.blacklist:
								if field in not_required:
									datum = CONFIG.not_required_label
								elif field=="compound_type":
									#Attempt to guess missing compound_types.
									datum=guess_type(model_fields["compound"])
								elif allow_unknowns:
									datum = CONFIG.unknown_label ###Take this value or no?
									data_is_valid = False
								elif field=="compound" and not datum:
									raise Exception("No compound specified!")
								elif field=="abbrev":
									#Allow absent abbreviations, but mark them the same as the datum.
									no_abbrev=True
								else:
									raise Exception("Value not allowed: \"{}\"".format(datum))
							else:
								#If the field is a quantity, check for units.
								if field[:-2]=="quantity":
									#Remove punctuation and whitespace if necessary.
									datum = str(datum).lower().translate(None, "\n?/,!@#$%^&*-+=_\\|")

									#Make sure parentheses do not contain numbers as well -- and if they do, remove them.
									try:
										paren_contents = datum[datum.index("(")+1:datum.index(")")]
										if paren_contents.isalpha():
											unit = paren_contents
										else:
											#Ignore the paren_contents if obviously not a unit.
											datum = datum[:datum.index("(")]+datum[datum.index(")")+1:]
											###Display error or auto-convert?
											###raise Exception("Invalid character found in unit: {}".format(paren_contents))
									except:
										unit = "" #Gather a unit from quantity if present.

									stripped_datum = ""
									for element in datum:
										if element in "1234567890. ": #Ignore spaces as well.
											stripped_datum += element
										else:
											old_datum = datum
											datum = float(stripped_datum)

											#If another element is reached, assume it is a unit.
											#	(if it isn't valid, raise an exception)
											if element == "g": unit = "g"
											elif element == "m": unit = "mL"
											elif element == "d": unit = "d"
											else: raise Exception("Unknown unit present: {}".format(old_datum[old_datum.index(element):]))
											break #Ignore anything beyond the unit.

									#If the unit was auto-added_quantity, add the unit to the correct field.
									corresponding_unit = "unit_{}".format(field[-1])
									if corresponding_unit in auto_added_fields:
										#If no unit was gathered, default to grams.
										if unit=="": unit = "g"
										#Apply the unit to the data entry.
										model_fields[corresponding_unit] = unit

								#Trim "note" fields that are over the range.
								elif field=="notes":
									if len(datum) > int(data_ranges["notes"][1]):
										datum = datum[:int(data_ranges["notes"][1])]
								#Translate any yes/no answer to "Yes"/"No"
								elif field in bool_fields:
									datum = datum.lower()
									if "y" in datum: datum="Yes"
									elif "n" in datum: datum="No"
									elif "" in datum: datum="No"
								#Make sure references are not repeated.
								elif field=="ref":
									try:
										assert(not datum in ref_set)
									except:
										raise Exception("Reference already in use!")


								elif field == "CAS_ID":
									datum = datum.replace("/","-").replace(" ","-").replace("_","-")
								elif field == "abbrev":
									try:
										assert datum not in abbrev_dict
										abbrev_dict[datum] = True #Remember the abbreviation is now being "used."
									except:
										raise Exception("Abbreviation already used!")
								elif field == "compound_type":
									###Gross? Why not use dict, Past Casey? --Future Casey
									print "yes"
									for option in edit_choices["typeChoices"]:
										datum = datum[0:2].lower()
										if option[0:len(datum)].lower() == datum:
											print "Found match!"
											datum = option
											break
								try:
									#Attempt to validate the data.
									if datum != CONFIG.unknown_label:
										assert(quick_validation(field, datum, model=model))

								except:
									data_is_valid = False
									raise Exception("Datum did not pass validation!")

								try:
									#Post-validation configuration.
									if model=="CompoundEntry":
										if field == "compound":
											try:
												chemspider_data = chemspipy.find_one(datum)
												#print "RAW:\t{}".format(datum)
												#print "ChemSpi:\t{}: {}\n".format(chemspider_data.commonname, chemspider_data.smiles)
											except:
												raise Exception("Unknown compound! Try another name?")
										if no_abbrev and field=="compound":
											print "No abbreviation found for \"{}\" but continuing!".format(datum)###
											model_fields["abbrev"] = datum
								except:
									raise Exception("Post-validation CONFIGuration failed!")

							model_fields[field] = datum
							#Continue to iterate through the row.
							i+=1
							j+=1
						except Exception as e:
							raise Exception([row_num, str(e)])

					#Add the new entry to the list of entries to batch add.
					if model=="Data":
						model_fields["is_valid"] = data_is_valid
						entry_list.append(new_Data_entry(u, **model_fields))
					elif model=="CompoundEntry":
						entry_list.append(new_CG_entry(lab_group, **model_fields))

					added_quantity += 1
				except Exception as e:
					if type(e.args[0])==list:
						error_log.append(e.args[0]+[field, datum])
					else:
						error_log.append([str(e)])
					error_quantity +=1
				row_num += 1

			#Update cursors if data was added.
			if added_quantity and model=="Data":
				Data.objects.bulk_create(entry_list)

				#Remove the cached version of the last page.
				old_data_size = get_cache(lab_group, "TOTALSIZE")
				control.clear_page_of(lab_group, old_data_size)

				#Clear any other caches.
				set_cache(lab_group, "REFS", None)

				#Add the new data to the cached size.
				set_cache(lab_group, "TOTALSIZE", old_data_size+len(entry_list))

			elif added_quantity and model=="CompoundEntry":
				CompoundEntry.objects.bulk_create(entry_list)
				set_cache(lab_group, "COMPOUNDGUIDE", None)
				set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", None)
			entry_list = None #Clear the entry_list out of memory.
		except Exception as e:
			error_log.append([str(e)])
	elif not u.is_authenticated():
		fatal_message = "Please log in to upload data."
		return HttpResponse(fatal_message)
	elif request.method=="POST" and not request.FILES:
		fatal_message = "No file uploaded."
	else:
		return render(request, 'upload_form.html')

	if error_quantity != 0:
		success_percent = "{:.1%}".format(float(added_quantity)/(added_quantity+error_quantity))
	else:
		success_percent = "100%"

	#Cache the error_log for easy access later.
	if error_log:
		#Store the error log for 4 hours.
		set_cache(u.get_profile().lab_group, "UPLOADERRORLOG", 14400)
	#Render the results template.
	return render(request, 'upload_results.html', {
		"fatal_message": fatal_message, #Includes data and data_indexes.
		"error_log": error_log,
		"added_quantity": added_quantity,
		"error_quantity": error_quantity,
		"success_percent": success_percent,
	})

def download_CSV(request): ###Need to fix.
	u = request.user
	#Make sure the user is authenticated before downloading data.
	if not u.is_authenticated():
		return HttpResponse("<p>Please log in to download data.</p>")

	if request.method=="POST":
		#Specify which CSV to download.
		model = request.POST["dataType"]
		data_filter = request.POST["downloadFilter"]

		#Generate a file name.
		lab_group = u.get_profile().lab_group
		date = datetime.datetime.now()
		file_name = "{}_{}_{}".format(lab_group.lab_title, model, data_filter).replace(" ", "_").lower()

		#Set up the HttpResponse to be a CSV file.
		CSV_file = HttpResponse(content_type="text/csv")
		CSV_file["Content-Disposition"] = "attachment; filename={}.csv".format(file_name)

		#Django HttpResponse objects can be handleded like files.
		writer = csv.writer(CSV_file)

		#Write the verbose headers to the CSV_file
		verbose_headers = get_model_field_names(verbose=True, model=model) ###NOT TESTED.
		writer.writerow(verbose_headers)

		#Write the actual entries to the CSV_file if the user is authenticated.
		headers = get_model_field_names(verbose=False, model=model)

		if model=="Data":
			collected_data = control.collect_all_data(lab_group)
			if data_filter=="good":
				CSV_data = collected_data.filter(is_valid=True)
			elif data_filter=="bad":
				CSV_data = collected_data.filter(is_valid=False)
			else:
				CSV_data = collected_data

		elif model=="CompoundEntry":
			###Better way? Generalize mass CG collection?
			CSV_data = CompoundEntry.objects.filter(lab_group=lab_group).order_by("compound")

		for entry in CSV_data:
			try:
				#Create and apply the actual row.
				row = [getattr(entry, field).encode("utf-8") for field in headers]
				writer.writerow(row)
			except Exception as e:
				print("ERROR getting '{}':{}".format(entry,e))###
		return CSV_file #ie, return HttpResponse(content_type="text/csv")
	else:
		return render(request, 'download_form.html')

def download_error_log(request): ###Nothing done yet... ;B
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
		verbose_headers = get_model_field_names(verbose=True)
		writer.writerow(verbose_headers)

		#Write the actual entries to the CSV_file if the user is authenticated.
		headers = get_model_field_names(verbose=False)
		lab_data = control.collect_all_data(u.get_profile().lab_group)

		for entry in lab_data:
			row = []
			try:
				#Apply the other columns.
				for field in headers:
					if field[-1].isdigit():
						pass
					else:
						row += [eval("entry.{}".format(field))]
				writer.writerow(row)
			except Exception as e:
				print(e)###
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
			session = control.get_page_info(request, int(num))

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
		lab_group = u.get_profile().lab_group
		name_pairs = collect_CG_name_pairs(lab_group, overwrite=False)
		return HttpResponse(json.dumps(name_pairs), mimetype="application/json")
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

		while (len(changesMade["edit"]) > 0):
			try:
				#An editPackage is [indexChanged, fieldChanged, newValue].
				editPackage = changesMade["edit"].pop()
				indexChanged = int(editPackage[0])-1 #Translate to 0-based Index
				fieldChanged = editPackage[1]
				newValue = editPackage[2] #Edits are validated client-side.
				dataChanged = lab_data[indexChanged]

				#Specific field-checks.
				if fieldChanged=="ref":
					#If the value is valid, change all of the duplicated reactions
					#	to the new reference.
					assert(not newValue in get_ref_set(lab_group))
					old_ref = getattr(dataChanged, fieldChanged)
					Data.objects.all().filter(lab_group=lab_group, duplicate_of=old_ref).update(duplicate_of=newValue)
					control.clear_all_page_caches(lab_group)
				else:
					control.clear_page_of(lab_group, indexChanged)

				#Make the edit in the database
				setattr(dataChanged, fieldChanged, newValue)
				#Data should only move closer to being "valid data."
				if not dataChanged.is_valid:
					revalidate_data(dataChanged, lab_group)
				dataChanged.user = u
				dataChanged.save()

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

				control.clear_page_of(lab_group, total_size)
				###Need to change size?
			except:
				pass
		while (len(changesMade["del"]) > 0):###SLOW
			try:
				indexChanged = changesMade["del"].pop()-1 #0-based Index
				lab_data[indexChanged].delete() #Django handles the deletion process.
				control.clear_page_of(lab_group, indexChanged)
			except:
				pass

	return HttpResponse("OK"); #Django requires an HttpResponse...

def get_full_datum(request, control=control):
	u = request.user
	try:
		if u.is_authenticated() and request.method=="POST":
			#Give the specific index requested.
			index_requested = json.loads(request.POST["indexRequested"])
			entry = control.collect_all_data(u.get_profile().lab_group)[index_requested]
			return render(request, "full_datum.html", {
				"entry":entry,
			})
		else: return HttpResponse("Please log in to see data details.")

	except Exception as e:
		print(e)
		return HttpResponse("<p>Data could not be loaded!</p>")

######################  User Auth ######################################
def change_password(request):
	code_sent = False
	if request.method == "POST":
		print(request.POST) ####
		if request.POST.get("emailCode"):
			print("YES\n")
			code_sent = True
		username = request.POST.get("email", "")
		code = request.POST.get("code", "")
		password = request.POST.get("newpassword", "")

	return render(request, "change_password_form.html", {
		"code_sent":code_sent
	})

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
	return HttpResponse("OK")

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

######################  Developer Functions  ###########################

######################  Error Messages  ################################
def display_404_error(request):
	response = render(request, '404_error.html')
	response.status_code = 404
	return response

def display_500_error(request):
	response = render(request, '500_error.html')
	response.status_code = 500
	return response
