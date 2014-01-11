from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.core.mail import send_mail
from django.contrib.auth.hashers import *
from django.shortcuts import render
from django.contrib import auth
from django.db.models import Q
from models import *
from forms import *
from validation import *

import json
import csv
import string
import datetime
import rdkit.Chem as Chem

#from svg_construction import *
#from construct_descriptor_table import *
from data_config import CONFIG

# # # # # # # # # # # # # # # # # # #
  # # # # # # # # Data and Page Helper Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
def get_lab_data(lab_group):
 return Data.objects.filter(lab_group=lab_group).order_by("creation_time")

def get_public_data():
 #Only show the public data that is_valid.
 return Data.objects.filter(public=True, is_valid=True).order_by("creation_time")

#Returns a specific datum if it is public or if it belongs to a lab_group.
def get_datum(lab_group, ref):
 query = Data.objects.filter(Q(ref=ref), Q(lab_group=lab_group) | Q(public=True))
 if not query.exists():
  raise Exception("Datum not found!")
 return query[0]
 
def get_lab_data_size(lab_group):
 size = get_cache(lab_group, "TOTALSIZE")
 if not size:
  size = get_lab_data(lab_group).count()
  set_cache(lab_group, "TOTALSIZE", size)
 return size

def calc_total_pages(data_size):
 total_pages = 1 + int((data_size-1)/CONFIG.data_per_page)
 return total_pages if total_pages > 1 else 1

def get_page_link_format(current, total):
 #Variable Setup:
 radius = CONFIG.current_page_radius
 data_per_page = CONFIG.data_per_page

 #Always display at least page one.
 page_set = {1}
 if total==1:
  return list(page_set)
 else:
  for i in xrange(current-radius, current+radius+1):
   if (1 < i < total): page_set.add(i)
  page_set.add(total) #Always show the last page.
 
  #Convert page_links to an ordered list.
  raw_page_links = list(page_set)
  raw_page_links.sort()

  page_links = []
  for i in raw_page_links:
   if page_links and i-page_links[-1]>1:
    page_links.append("...")#Add an ellipsis between any pages that have a gap between them.
    page_links.append(i)
   else:
    page_links.append(i) 
  return page_links

def get_pagified_data(page, lab_group=None, data=None):
 #Variable Setup:
 data_per_page = CONFIG.data_per_page

 #Get the Lab's data if data is not specified. 
 if not data:
  try:
   data = get_lab_data(lab_group)
  except:
   raise Exception("No data nor lab_group was specified.")
 total_pages = calc_total_pages(data.count())

 #Check that the page request is valid.
 if not (0<page<=total_pages):
  raise Exception("Page requested is outside of page range.")

 #Return the data that would be contained on the requested page.
 page_data = data[(page-1)*data_per_page:page*data_per_page]
 return page_data

#Get data before/after a specific date (ignoring time).
def get_date_filtered_data(lab_group, raw_date, direction="after", lab_data=None):
 #Convert the date input into a usable string. (Date must be given as MM-DD-YY.)
 date = str(datetime.datetime.strptime(raw_date, "%m-%d-%y"))

 #Only get the data that belongs to a specific lab_group.
 if lab_data:
  lab_data = lab_data.filter(lab_group=lab_group)
 else:
  lab_data = get_lab_data(lab_group)

 #Get the reactions before/after a specific date.
 if direction.lower() == "after":
  filtered_data = lab_data.filter(creation_time__gte=date_string)
 else:
  filtered_data = lab_data.filter(creation_time__lte=date_string)

 return filtered_data
 
#Returns the info that belongs on a specific page.
def get_page_info(request, page = None, data=None):
 try:
  #Gather necessary information from the user's session:
   #Get the user's current_page (or 1 if it is unknown.
  if not page:
    page = int(request.COOKIES.get("current_page", 1))
  
  if not data:
   u = request.user
   if u.is_authenticated():
    data = get_lab_data(u.get_profile().lab_group) ###TODO: Union this with public_data for lab_group?
   else:
    data = get_public_data()
  
  total_data_size = data.count()
  total_pages = calc_total_pages(total_data_size)

  #Make sure the page is a valid page. If not, go to the last page.
  if not (0 < page <= total_pages):
   page = total_pages

  #Pack up the session info:
  session = {
   "page_data": get_pagified_data(page, data=data),
   "total_data_size": total_data_size,
   "total_pages": total_pages,
   "page_links": get_page_link_format(page, total_pages),
   "current_page": page,
  }
  return session
 except Exception as e:
  raise Exception("Data could not be retrieved for page {}:\n--{}.".format(page, e))

def repackage_page_session(session):
 data = session["page_data"]
 total_pages = session["total_pages"]
 page_links = session["page_links"]
 current_page = session["current_page"]
 total_data_size = session["total_data_size"]

 #Show the overall index of each datum.
 start_index = (current_page-1)*CONFIG.data_per_page + 1
 end_index = (current_page)*CONFIG.data_per_page + 1

 #Prepare packages.
 data_package = zip(data, range(start_index, end_index))
 page_package = {
  "current_page":current_page,
  "total_pages":total_pages,
  "data_per_page":CONFIG.data_per_page,
  "page_links":page_links,
  }
 return data_package, page_package, total_data_size

# # # # # # # # # # # # # # # # # # #
  # # # # # # # # Caching Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
def clear_page_cache(lab_group, page):
 set_cache(lab_group, "PAGEDATA|{}".format(page), None)

def clear_page_cache_of_index(lab_group, indexChanged):
 page = (indexChanged/CONFIG.data_per_page) + 1
 clear_page_cache(lab_group, page)

def clear_all_page_caches(lab_group, skip_data_check=False):
 if skip_data_check:
  total_pages = get_cache(lab_group, "TOTALPAGES")
 else:
  total_pages = calc_total_pages(get_lab_data(lab_group).count())
 for i in xrange(1, total_pages+1):
  clear_page_cache(lab_group, i)

# # # # # # # # # # # # # # # # # # #
  # # # # # # # # View Helper Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
def get_template_form(entry, model):
 result = {}
 if model=="Recommendation":
  result["reaction"] = [[getattr(entry, "reactant_{}".format(i)), getattr(entry, "quantity_{}".format(i)), getattr(entry, "unit_{}".format(i))] for i in CONFIG.reactant_range()]
  fields = get_model_field_names(model="Recommendation", unique_only=True)
  verbose_fields = get_model_field_names(verbose=True, model="Recommendation", unique_only=True)
  result["info"] = [[i, getattr(entry,j)] for (i,j) in zip(verbose_fields, fields)] 
 else:
  raise Exception("Model type not found!") 
 return result

#Given a user, change their password and email them the new password.
def randomize_password(user):
 new_pass = get_random_code(15) #Generate a random password for the user.
 user.password = make_password(new_pass) #Hash the password. 
 user.save()
 email_body = "Hello {},\n\n According to our records, you just requested a password change. We have changed your account information as follows:\nUsername: {}\nPassword: {}".format(user.first_name, user.username, new_pass)
 send_mail("Dark Reactions: Password Change Request", email_body, settings.EMAIL_HOST_USER, [user.email], fail_silently=False)

# # # # # # # # # # # # # # # # # # #
  # # # # # # # # View Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
def database(request, global_data=False):
 #Organize the session information.
 session = get_page_info(request)
 data_package, page_package, total_data_size = repackage_page_session(session)
 return render(request, 'global_page.html', {
  "data_on_page": data_package, #Includes data and data_indexes.
  "page_package": page_package,
  "total_data_size": total_data_size,
  "template": "database",
 })

def data_transmit(request):
 try:
  try:
   #Variable Setup
   u = request.user
   body = json.loads(request.POST["body"], "utf-8")
   query_list = body.get("currentQuery")
   page = body.get("page")

   if query_list:
    data = filter_data(u.get_profile().lab_group, query_list)
   else:
    data = get_lab_data(u.get_profile().lab_group)
   session = get_page_info(request, page=int(page), data=data)
  except Exception as e:
   print e
   return HttpResponse("Page could not be loaded.")

  #Organize the data for the template.
  data_package, page_package, total_data_size = repackage_page_session(session)
  return render(request, 'database.html', {
   "data_on_page": data_package, #Includes data indexes
   "page_package": page_package, #Includes page links
   "total_data_size": total_data_size,
  })
 except Exception as e:
  print e
  return HttpResponse("Page \"{}\" could not be loaded".format(page))

def recommend(request): ###TODO: ADD TEMPLATE BITS, CASEY! 
 #Get user data if it exists.
 u = request.user
 if u.is_authenticated():
  fatal_message = ""
  recommendation_query = get_recommendations_by_date(u.get_profile().lab_group)
  recommendations = [get_template_form(i, "Recommendation") for i in recommendation_query] 
 else:
  fatal_message = "Please log in to view recommendations."
  recommendations = None

 return render(request, 'recommend_global.html', {
  "recommendations": recommendations,
  "fatal_message": fatal_message,
 })


###TODO: Add this -- right now it is just a copy. 
def saved(request): ###TODO: ADD TEMPLATE BITS, CASEY! 
 #Get user data if it exists.
 u = request.user
 if u.is_authenticated():
  fatal_message = ""
  recommendation_query = get_recommendations_by_date(u.get_profile().lab_group)
  recommendations = [get_template_form(i, "Recommendation") for i in recommendation_query] 
 else:
  fatal_message = "Please log in to view recommendations."
  recommendations = None

 return render(request, 'recommend_global.html', {
  "recommendations": recommendations,
  "fatal_message": fatal_message,
 })


# # # # # # # # # # # # # # # # # # #
  # # # # # # # # Sub-view Functions (eg, Javascript response views) # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
def get_user_license_agreement(request):
 u = request.user
 #Indicate whether the user needs to agree to updated terms or sign the terms initially.
 if u.is_authenticated():
  if u.get_profile().license_agreement_date:
   if u.get_profile().license_agreement_date < CONFIG.current_license_date:
    license_changed = True
   else:
    return HttpResponse("Our records indicate that you already agreed to the terms and conditions!")
  else:
   license_changed = False
   
  return render(request, 'user_license_form.html', {
   "license_changed": license_changed,
   "license_file": CONFIG.current_license_file,
   "license_date": CONFIG.current_license_date.split(" ")[0], #Show the modification date but ignore the time.
  })
 else:
  return HttpResponse("Please create a user to accept the license agreement.")

def update_user_license_agreement(request):
 u = request.user 
 if request.method=="POST":
  try:
   u.get_profile().license_agreement_date = str(datetime.datetime.now())
   u.get_profile().save()  
   return HttpResponse("Thank you for joining the Dark Reactions Project." +
   "<button class=\"genericButton form_button clearScreenButton\">Explore</button>")
  except:
   return HttpResponse("Your request could not be completed. Please try again.")
 else:
  return HttpResponse("Please click the \"I Agree\" to accept the Terms and Conditions.")

def save_recommmendation(request):
 u = request.user
 if request.method=="POST":
  try:
   lab_group = u.get_profile().lab_group
   rec_info = json.loads(request.POST.get("rec_info"))
   if not rec["date"]:
    #If no date is specified, assume the latest data is what is being changed.
    rec["date"] = get_recommendations_by_date(lab_group).get().date  
   rec = Recommendation.objects.filter(date=date).get()
   
  except:
   return HttpResponse("Your request could not be completed. Please try again.")
 else:
  return HttpResponse("Please log in to save recommendations")
  
####################################################
####################################################
####################################################
####################################################
############   REWRITTEN THUS FAR ##################
####################################################
####################################################
####################################################
####################################################
def revalidate_all_data(lab_group, invalid_only = True, return_errors=False):
 data_to_validate = Data.objects.filter(lab_group=lab_group)

 if invalid_only:
  data_to_validate = data_to_validate.filter(is_valid=False)

 print "Found {} to validate.".format(data_to_validate.count())
 error_log = {"CG":[], "quantity": [], "refs":[], "pH":[], "temp":[], "time": []}
 j = 0

 if data_to_validate.exists():
  for data in data_to_validate:
   #Display how many were analyzed so far.
   if j%500==0: print "Played with {}...".format(j)

   errors = revalidate_data(data, lab_group, batch=True)
   for i in CONFIG.reactant_range():
    if "reactant_{}".format(i) in errors:
     error_log["CG"].append("{}".format(getattr(data, "reactant_{}".format(i))))
   for i in CONFIG.reactant_range():
    if "quantity_{}".format(i) in errors:
     error_log["quantity"].append(getattr(data, "ref"))
   if "temp" in errors:
    error_log["temp"].append(getattr(data, "ref"))
   if "time" in errors:
    error_log["time"].append(getattr(data, "ref"))
   if "pH" in errors:
    error_log["pH"].append(getattr(data, "ref"))

  #Clear the page caches
  clear_all_page_caches(lab_group)
  if return_errors:
   return error_log


######################  Core Views  ####################################

###TODO: Modify to get recommendations from DB
def get_reactants(lst): ###TODO: GENERALIZE
 reactant_info = [[lst[i-3],lst[i-2],lst[i-1]] for i in range(3,len(lst),3)]
 return reactant_info
  
def get_info(lst): ###TODO: GENERALIZE
 #Get from Temperature on -- excluding Notes.
 verbose_headers = get_model_field_names(verbose=True, model="Data")[16:]
 info = [[i, j] for (i,j) in zip(verbose_headers, lst)]
 return info 

def get_recommendations_by_date(lab_group, date = "recent"):
 if date=="recent":
  #Get the most recent version of the model.
  try:
   version = Model_Version.objects.filter(lab_group=lab_group, model_type="Recommendation").order_by("date")[0]
   date = version.date
  except Exception as e:
   raise Exception("Could not find any version of the model: {}".format(e))

 #Get the data associated with a specific date. ###TODO: Probably should be "creation_time"
 try:
  recommendations = Recommendation.objects.filter(lab_group=lab_group, date=date).order_by("score")
 except Exception as e:
  raise Exception("Could not find any version of the model: {}".format(e))

 return recommendations


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
    print "Generating SVG..."
    svg = generate_svg(u.get_profile().lab_group)
    print "Took {} seconds overall.".format(time.time()-start_time)###

    set_cache(lab_group, "TESTSVG", svg)
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
  #Rules:
  # 1.) A Lab can search only its data and public data.
  # 2.) Data that is returned can be sent to editing functions.

  #Verbose JSON Formats:
  # request.body ===
  #  query_list --> {[{u"field":u"FIELD", u"value":u"VALUE"}, ...]}

def filter_data(lab_group, query_list):
 #Variable Setup:
 lab_data = get_lab_data(lab_group)
 filters = ""

 #Collect all the valid search options
 non_reactant_fields = get_model_field_names(unique_only=True)
 legal_fields = set(non_reactant_fields+["reactant","quantity","unit","public","is_valid", "user"])

 #Check the query_list input before performing any database requests.
 try:
  for query in query_list:
   assert query.get(u"field") in legal_fields
   assert query.get(u"value")
   assert not "\"" in query.get(u"value")

 except:
  raise Exception("One or more inputs is illegal.")
 
 try:
  for query in query_list:
   field = query.get(u"field")
   value = query.get(u"value")
 
   if field in list_fields:
    #Check all the reactant/quantity/unit fields.
    Q_obj = ''.join(["Q({}_{}__icontains=\"{}\")|".format(field, i, value) for i in CONFIG.reactant_range()])[:-1]
    filters += ".filter({})".format(Q_obj)
   elif field=="atoms":
    atom_list = value.split(" ")
    if len(atom_list)>1:
     search_bool = atom_list.pop(-2) #Take the "and" or "or" from the list.
     op = "," if search_bool == "and" else "|" #Assign the correct Q operator.
     #Add the atoms  to a Q filter.
     Q_obj = ''.join(["Q(atoms__contains=\"{}\"){}".format(atom, op) for atom in atom_list])[:-1]
    else:
     Q_obj = "Q(atoms__contains=\"{}\")".format(atom_list[0])
    filters += ".filter({})".format(Q_obj)
   else:
    #Translate Boolean inputs into Boolean values.
    if field in bool_fields:
     value = True if value.lower()[0] in "1tyc" else False
     filters += ".filter({}={})".format(field, value)
    else:
     filters += ".filter({}__icontains=\"{}\")".format(field, value)
  data = eval("lab_data"+filters).order_by("creation_time")
  return data

 except Exception as e:
  print e
  raise Exception("Woops! A problem has occurred...")

def search(request):
 u = request.user

 #Collect the fields that will be displayed in the Search "Fields" tab.
 search_fields = get_model_field_names(both=True, unique_only=True)
 search_fields = [
  {"raw":"reactant", "verbose":"Reactant"},
  {"raw":"quantity", "verbose":"Quantity"},
  {"raw":"unit", "verbose":"Unit"},
  {"raw":"is_valid", "verbose":"Is Valid"},
  {"raw":"user", "verbose":"User"},
  {"raw":"public", "verbose":"Public"}] + search_fields

 if u.is_authenticated() and request.method=="POST":
  body = json.loads(request.POST["body"], "utf-8")
  query_list = body.get("currentQuery")
  
  try:
   lab_group = u.get_profile().lab_group
   #If no query_list is given, just return the lab_data.
   if query_list:
    data = filter_data(lab_group, query_list)
   else:
    data = get_lab_data(lab_group)
  except Exception as e:
   return HttpResponse(e)  

  #Return the filtered data to the user.
  return data_transmit(request)

 else:
  return render(request, 'search_global.html', {
  "search_fields": search_fields,
  })

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

def compound_guide_entry(request):
 u = request.user
 if u.is_authenticated():
  lab_group = u.get_profile().lab_group
  if request.method == 'POST':
   entry_info = json.loads(request.body, "utf-8")
   query = CompoundEntry.objects.filter(lab_group=lab_group, compound=entry_info["compound"])
   if query.exists():
    return render(request, 'compound_guide_row.html', {
     "entry": query[0]
    })
   else:
    return HttpResponse("No CG found. Please refresh page.")
  else:
   return HttpResponse("Please use the compound guide interface.")
 else:
  return HttpResponse("Please log in to access the compound guide!")

def edit_CG_entry(request): ###Edits?
 u = request.user
 if request.method == 'POST' and u.is_authenticated():
  changesMade = json.loads(request.body, "utf-8")

  #Get the Lab_Group data to allow direct manipulation.
  lab_group = u.get_profile().lab_group
  CG_data = collect_CG_entries(lab_group)

  if changesMade["type"]=="del":
   for i in changesMade["data"]:
    try:
     if i["abbrev"]=="": #If the compound has no abbrev.
      CompoundEntry.objects.filter(lab_group=lab_group, abbrev=i["compound"], compound=i["compound"])[0].delete()
     else:
      CompoundEntry.objects.filter(lab_group=lab_group, abbrev=i["abbrev"], compound=i["compound"])[0].delete()
    except Exception as e:
     print("Could not delete!\n {}".format(e))
  elif changesMade["type"]=="edit":
   try:
    #Variable Setup
    field = changesMade["field"]
    new_val  = changesMade["newVal"]
    old_val  = changesMade["oldVal"]
    compound = changesMade["compound"]

    assert(compound and field)

    #Gather all relevant data (and check for duplicates).
    try:
     changed_entry = CompoundEntry.objects.get(lab_group=lab_group, compound=compound)
    except:
     return HttpResponse("Please delete duplicate entry!")

    #Make sure the compound isn't already being used.
    if field=="compound":
     possible_entry = CompoundEntry.objects.filter(lab_group=lab_group, compound=new_val)
     if possible_entry.exists() and possible_entry[0].compound!=compound:
      return HttpResponse("Already used!")
    elif field=="abbrev":
     possible_entry = CompoundEntry.objects.filter(lab_group=lab_group, abbrev=new_val)
     if possible_entry.exists() and possible_entry[0].compound!=compound:
      return HttpResponse("Already used!")
    elif field=="CAS_ID" and new_val: #Empty CAS_IDs don't require queries.
     possible_entry = CompoundEntry.objects.filter(lab_group=lab_group, CAS_ID=new_val)
     if possible_entry.exists() and possible_entry[0].compound!=compound:
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

    #Lookup fresh data from ChemSpider and RDKit
    update_compound(changed_entry)


    #Make any edits to the Data if needed.
    if field=="abbrev":
     #Don't overwrite ALL empty entries -- only those related to the compound.
     if not old_val:
      old_val = compound

     for i in CONFIG.reactant_range():
      exec("old_data = Data.objects.filter(lab_group=lab_group, reactant_{}=old_val)".format(i))
      exec("old_data.update(reactant_{}=new_val)".format(i))

    #If no inorganic was found, ask the user for its identity.
    if changed_entry.compound_type=="Inorg":
     return HttpResponse("Not found!")

    ###TODO: Clear the cached CG
   except Exception as e:
    return HttpResponse("Need more info!")

  #Clear the cached CG entries.
  set_cache(lab_group, "COMPOUNDGUIDE", None)
  set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", None)

 return HttpResponse("0")

 def add_inorg_info(request):
	 pass#Add "smiles" and "mw" info.

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
def data_form(request, copy_ref=None): #If no data is entered, stay on the current page.
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
   old_data_size = get_lab_data_size(lab_group)
   clear_page_cache_of_index(lab_group, old_data_size)

   #Refresh the TOTALSIZE cache.
   set_cache(lab_group, "TOTALSIZE", old_data_size + 1)
   success = True #Used to display the ribbonMessage.
 else:
  #Submit a blank form or an auto-filled form if a ref is supplied.
  if copy_ref==None:
   form = DataEntryForm(
    initial={"leak":"No"}
   )
  elif u.is_authenticated():
   datum = get_datum(u.get_profile().lab_group, copy_ref.replace("+"," "))
   initial_fields = {field:getattr(datum, field) for field in get_model_field_names(model="Data")}
   form = DataEntryForm(
    initial=initial_fields
   )
  else:
   return HttpResponse("You do not have access to copy this datum.")
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
  elif "dupl" in stripped_heading:
   related_field = "duplicate_of"
  elif "recom" in stripped_heading:
   related_field = "recommended"
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
      # and that all fields are valid.
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
       # no data is present for the generated column.
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
           # (if it isn't valid, raise an exception)
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
         for option in edit_choices["typeChoices"]:
          datum = datum[0:2].lower()
          if option[0:len(datum)].lower() == datum:
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
            ###The query is sent here for validation, AND later for data; super inefficient. --Casey (TODO)
            chemspider_data = chemspipy.find_one(datum)
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
      model_fields["image_url"], model_fields["smiles"], model_fields["mw"] = chemspider_lookup(model_fields)
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
    clear_page_cache_of_index(lab_group, old_data_size)

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

  #Specify which data filter was used.
  data_filter = request.POST["downloadFilter"]
  if data_filter == "simple":
   data_filter = ""

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


  #Write the actual entries to the CSV_file if the user is authenticated.
  headers = get_model_field_names(verbose=False, model=model)

  if model=="Data":
   #Throw the "Ref" in front of the data if applicable.
   verbose_headers.remove("Reference")
   headers.remove("ref")
   verbose_headers.insert(0, "Reference")
   headers.insert(0, "ref")

   #Write the header row
   writer.writerow(verbose_headers)

   #Gather the data and filter it depending on what the user wants.
   collected_data = get_lab_data(lab_group)

   if data_filter=="good":
    CSV_data = collected_data.filter(is_valid=True)
   elif data_filter=="bad":
    CSV_data = collected_data.filter(is_valid=False)
   else:
    CSV_data = collected_data

  elif model=="CompoundEntry":
   if data_filter=="complex":
    verbose_headers.append("SMILES")
    headers.append("smiles")

   #Write the header row
   writer.writerow(verbose_headers)

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
  lab_data = get_lab_data(u.get_profile().lab_group)

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
  #Rules:
  # 1.) A Lab can only delete data it owns.
  # 2.) Users can only modify their own Lab's data.

  #Verbose JSON Formats:
  # request.body ===
  #  delete_Data --> {[originalRef_1, ..., originalRef_N]}
  #  change_Data --> {originalRef:"X", fieldChanged:"X", newValue:"X"}

def delete_Data(request):
 u = request.user
 if request.method == "POST" and u.is_authenticated():
  #Variable Setup
  lab_group = u.get_profile().lab_group
  deleteList = json.loads(request.body, "utf-8")
  lab_data = get_lab_data(lab_group)

  #Find and delete data entries in a User's Lab. 
  for ref in deleteList:
   try:
    lab_data.filter(ref=ref).first().delete()
   except:
    HttpResponse("One or more selected data not found.")

  #Finally, return a success code.
  return HttpResponse(0);
 return HttpReponse("Please log in to delete data.") 

def change_Data(request):
 #Fields that may be changed via this script.
 whitelist = set(get_model_field_names())

 u = request.user
 if request.method == "POST" and u.is_authenticated():
  #Variable Setup
  lab_group = u.get_profile().lab_group
  editLog = json.loads(request.body, "utf-8")
  lab_data = get_lab_data(lab_group)
  
  #Get the Datum for the lab.
  ref = editLog["ref"]
  fieldChanged = editLog["field"]
  newValue = editLog["newValue"]
  datum = lab_data.filter(ref=ref).first() 

  #Check that the field being edited actually exists (to prevent sabotage).
  if not fieldChanged in whitelist:
   return HttpResponse("Cannot modify {}".format(fieldChanged))

  #Check that the edit doesn't invalidate the Datum in any way.   
  try: 
   oldValue = getattr(datum, fieldChanged)
   setattr(datum, fieldChanged, newValue)
   dirty_data = model_to_dict(datum, fields=get_model_field_names())
   if fieldChanged=="ref":
    clean_data, errors = full_validation(dirty_data, lab_group)
   else:
    clean_data, errors = full_validation(dirty_data, lab_group, revalidating=True)
   
   #Send back an error if it exists.
   if errors:
    return HttpResponse("Error: {}".format(errors[errors.keys()[0]]))
 
   #Get the parsed value after cleaning.
   setattr(datum, fieldChanged, clean_data[fieldChanged])
 
   #Make the edit in the database.
   datum.user = u
   datum.save()
 
   if fieldChanged=="ref":
    #Update the "ref" in any Data of which it is a duplicate.
    lab_data.filter(duplicate_of=oldValue).update(duplicate_of=newValue)
 
   return HttpResponse(0)

  except:
   return HttpResponse("Edit unsuccessful...")

 return HttpResponse("Please log in to modify data.")

######################  User Auth ######################################
def change_password(request):
 error=False
 if request.method == "POST":
  try:
   email = request.POST.get("email")
   username = request.POST.get("username")
   last_name = request.POST.get("lastName")
   user = User.objects.filter(Q(email=email)|Q(username=username), Q(last_name=last_name))[0]
   #Change the user's password and send them an email.
   randomize_password(user)
   return HttpResponse("A new password has been emailed to you.")
  except:
   #If no user is found given the credentials, tell the user.
   error=True
 return render(request, "change_password_form.html", {
  "error":error
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
 return render(request, "registration_cell.html", {})

def user_update(request):
 u = request.user
 if request.method == "POST":
  form = UserForm(request.POST, instance=u)
  if form.is_valid():
   form.save()
   return HttpResponse("Update Successful!")
 else:
  form = UserForm(instance=u)
 return render(request, "user_update_form.html", {
  "form": form,
 })

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
  form = LabForm(data = request.POST)
  if form.is_valid():
   lab_group = form.save()
   access_code = lab_group.access_code
   
   #Send a "confirmation" email to the new lab email.
   email_body = "Thank you for joining the Dark Reactions Project!\n\nPlease continue by creating a \"user\" for your lab. Simply...\n\t1.) Record the \"access code\" for your lab: {}\n\t2.) Click \"Register\" and create a user using the access code above.\n\t3.) Start uploading data!\n\nWe wish you all the best,\nThe Dark Reactions Project Team".format(lab_group.access_code)

   send_mail("Dark Reactions: Lab Registration Successful", email_body, settings.EMAIL_HOST_USER, [lab_group.lab_email], fail_silently=False)

   #Send the DRP Admins an email about the new Lab Registration.
   #TODO:Remove this when we scale to unmanageable quantities of labs.
   alert_about_new_lab(lab_group)
 
   return HttpResponse("Registration Successful! Please check your email.")
 else:
  form = LabForm()
 return render(request, "lab_registration_form.html", {
  "form": form,
 })

######################  Developer Functions  ###########################
def alert_about_new_lab(lab_group):
 email_body = "A new Lab Group has been registered:"
 email_body += "\n{}\n{}\n{}".format(lab_group.lab_title, lab_group.lab_address, lab_group.lab_email) 
 send_mail("Dark Reactions: New Lab Group Registered", email_body, settings.EMAIL_HOST_USER, [settings.EMAIL_HOST_USER], fail_silently=False)

######################  Error Messages  ################################
def display_404_error(request):
 response = render(request, '404_error.html')
 response.status_code = 404
 return response

def display_500_error(request):
 response = render(request, '500_error.html')
 response.status_code = 500
 return response
