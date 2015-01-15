from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
from django.contrib.auth.hashers import *
from django.shortcuts import render

from database_construction import *
from forms import *
from validation import *

import json
import csv
import string

from data_config import CONFIG

# # # # # # # # # # # # # # # # # # #
  # # # # # # # # Data and Page Helper Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
def calc_total_pages(data_size):
 total_pages = 1 + int((data_size-1)/CONFIG.data_per_page)
 return total_pages if total_pages > 1 else 1


def get_pagified_data(page, lab_group=None, data=None):
 from DRP.models import get_lab_Data

 #Variable Setup:
 data_per_page = CONFIG.data_per_page

 #Get the Lab's data if data is not specified.
 if not data:
  try:
   data = get_lab_Data(lab_group)
  except:
   raise Exception("No data nor lab_group was specified.")
 total_pages = calc_total_pages(data.count())

 #Check that the page request is valid.
 if not (0<page<=total_pages):
  raise Exception("Page requested is outside of page range.")

 #Return the data that would be contained on the requested page.
 page_data = data[(page-1)*data_per_page:page*data_per_page]
 return page_data

#Returns the info that belongs on a specific page.
def get_page_info(request, page = None, data=None):
 from DRP.models import get_lab_Data
 from DRP.pagifier import get_page_link_format

 try:
  #Gather necessary information from the user's session:
   #Get the user's current_page (or 1 if it is unknown.
  if not page:
    page = int(request.COOKIES.get("current_page", 1))

  if not data:
   u = request.user
   if u.is_authenticated():
    data = get_lab_Data(u.get_profile().lab_group) ###TODO: Union this with public_data for lab_group?
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
 from DRP.models import get_lab_Data

 if skip_data_check:
  total_pages = get_cache(lab_group, "TOTALPAGES")
 else:
  total_pages = calc_total_pages(get_lab_Data(lab_group).count())
 for i in xrange(1, total_pages+1):
  clear_page_cache(lab_group, i)


# # # # # # # # # # # # # # # # # # #
  # # # # # # # # View Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
@login_required
@require_http_methods(["GET"])
def database(request, page_request=None):

  # If a valid page number was given, use it.
  try:
    page = int(page_request)
  except:
    page = None

  # Organize the session information.
  session = get_page_info(request, page=page)
  data_package, page_package, total_data_size = repackage_page_session(session)

  # Return a package of page information and data.
  return render(request, 'global_page.html', {
    "data_on_page": data_package, #Includes data and data_indexes.
    "page_package": page_package,
    "total_data_size": total_data_size,
    "template": "database",
  })

@login_required
@require_http_methods(["POST"])
def data_transmit(request):
 from DRP.models import get_lab_Data
 from DRP.retrievalFunctions import filter_data

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
    data = get_lab_Data(u.get_profile().lab_group)

   if data.count():
    session = get_page_info(request, page=int(page), data=data)
   else:
    return HttpResponse("No data found!")

  except Exception as e:
   print e
   return HttpResponse("Page could not be loaded.")

  #Organize the data for the template.
  data_package, page_package, total_data_size = repackage_page_session(session)
  return render(request, 'database.html', {
   "data_on_page": data_package, #Includes data indexes
   "page_package": page_package, #Includes page links
   "total_data_size": total_data_size,
   "template":"database",
  })
 except Exception as e:
  print e
  return HttpResponse("Page \"{}\" could not be loaded".format(page))

@login_required
def recommend(request, page_request=None):
  from DRP.pagifier import get_page_link_format

  # If a valid page number was given, use it.
  try:
    page = int(page_request)
  except:
    page = None

  #Get user data if it exists.
  u = request.user
  recs_per_page = 15

  recommendations = None
  fatal_message = ""
  try:
    recommendations = get_recommendations_by_date(u.get_profile().lab_group)
    total_data_size = recommendations.count()
    total_pages = total_data_size/recs_per_page

    if not request.GET.get("show_hidden"):
      recommendations = recommendations.filter(hidden=False)

    if not page:
      page = request.GET.get("rec_page", 1)

    recommendations = recommendations[(page-1)*recs_per_page:(page)*recs_per_page]

  except Exception as e:
    print e
    fatal_message = "No recommendations available."
    page = 1
    recs_per_page = 1
    page_links = []
    total_pages = 1

  return render(request, 'global_page.html', {
  "template":"recommendations",
  "recommendations": recommendations,
  "fatal_message": fatal_message,
  "total_data_size":total_data_size,
  "counter_offset":recs_per_page*(page-1),
  "page_package": {
                    "data_per_page":recs_per_page,
                    "current_page":page,
                    "page_links":get_page_link_format(page, total_pages),
                    }
 })


@login_required
@require_http_methods(["POST"])
def recommendation_transmit(request, seeded=False):
 from DRP.models import Recommendation

 try:
  #Variable Setup
  u = request.user
  body = json.loads(request.POST["body"], "utf-8")
  query_list = body.get("currentQuery")
  page = int(body.get("page", 1))

  #Get either the Seed Recs or the general Recommendations.
  if seeded:
    query_list.append({
	"field":"seeded",
	"value":"True",
	"match":"exact"})

  if query_list:
   recs = filter_recommendations(u.get_profile().lab_group, query_list)
  else:
   recs = get_recommendations_by_date(u.get_profile().lab_group)

  fatal_message = "" if recs.exists() else "No recommendations found."
  total_data_size = recs.count()


  recs_per_page = 10
  recs = recs[(page-1)*recs_per_page:(page)*recs_per_page]
  total_pages = total_data_size/recs_per_page

  template="seed_recommendations_entries.html" if seeded else "recommendations.html"
  return render(request, template, {
   "recommendations": recs,
   "fatal_message": fatal_message,
   "total_data_size":total_data_size,
   "counter_offset":recs_per_page*(page-1),
   "page_package": {
                    "data_per_page":recs_per_page,
                    "current_page":page,
                    "page_links":get_page_link_format(page, total_pages),
                    }
  })
 except Exception as e:
  print e
  return HttpResponse("Recommendations could not be loaded!")

@login_required
@require_http_methods(["POST"])
def edit_recommendation(request, action):
 try:
  u = request.user
  pid = request.POST["pid"]
  recommendations = get_recommendations(u.get_profile().lab_group)

  #Filter the recommendation of interest.
  rec = recommendations.get(id=pid)

  if action=="save":
   rec.saved = True
  elif action=="unsave":
   rec.saved = False
  elif action=="show":
   rec.hidden = False
  elif action=="hide":
   rec.hidden = True
  elif action=="sense":
   rec.nonsense = False
  elif action=="nonsense":
   rec.nonsense = True
  else:
   raise Exception("Error: Illegal action specified.")

  rec.user = u
  rec.save()

  return HttpResponse(0)
 except Exception as e:
  print e
  return HttpResponse(1)

@login_required
@require_http_methods(["GET"])
def saved(request):
 #Variable Setup
 u = request.user
 recommendations = None
 user_list = None

 fatal_message = ""
 try:
  #Get the recommendations that are saved for the lab.
  lab_group = u.get_profile().lab_group

  recommendations = get_recommendations(lab_group).filter(saved=True)
  assert recommendations.count()

  #Get the lab users for the select field.
  user_list = User.objects.filter(profile__lab_group=lab_group)

  #Get users
 except:
  fatal_message = "No saved recommendations available."

 return render(request, 'global_page.html', {
  "template":"saved",
  "recommendations": recommendations,
  "fatal_message": fatal_message,
  "users":user_list,
 })

@login_required
@require_http_methods(["GET"])
def rank(request):
 #Variable Setup:
 u = request.user
 unranked_rxn = get_random_unranked_reaction_or_none()

 if not unranked_rxn:
  fatal_message = "No unranked reactions available!"
 else :
  fatal_message = ""

 return render(request, 'global_page.html', {
  "template":"unranked_reaction",
  "unranked_rxn": unranked_rxn,
  "fatal_message": fatal_message,
 })


# # # # # # # # # # # # # # # # # # #
   # # # # # Sub-view Functions (eg, Javascript response views) # # # # # # # #
# # # # # # # # # # # # # # # # # # #
#Change the assigned_user of a recommendation.
@login_required
@require_http_methods(["POST"])
def assign_user_to_rec(request):
 u = request.user
 try:
  lab_group = u.get_profile().lab_group
  #Get PIDs from the request.
  rec_pid = request.POST["rec_pid"]
  user_pid = request.POST["user_pid"]

  #Get the recommendation to change.
  rec = get_recommendations(lab_group).get(id=rec_pid)

  #Get the newly assigned user or None.
  if user_pid:
   new_user = User.objects.get(id=user_pid)
  else:
   new_user = None

  #Assign the change in the database.
  rec.assigned_user = new_user
  rec.user = u
  rec.save()
  return HttpResponse(0)
 except Exception as e:
  print e
  return HttpResponse(1)

@login_required
@require_http_methods(["POST"])
def send_and_receive_rank(request):
 u = request.user
 #Get PIDs from the request.
 pid = request.POST["pid"]
 new_order = request.POST["newOrder"]

 #Get the recommendation to change.
 rxnlist = RankedReactionList.objects.get(id=pid)

 #Assign the change in the database.
 rxnlist.ranked_list = new_order
 rxnlist.ranker = u
 rxnlist.save()

 #Get a new (un)RankedReaction
 unranked_rxn = get_random_unranked_reaction_or_none()
 if not unranked_rxn:
  return HttpResponse("All unranked reactions now ranked!")

 return render(request, 'unranked_reaction_list.html', {
  "unranked_rxn": unranked_rxn,
 })

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
   return HttpResponse("<p>Your request could not be completed. Please try again.</p>")
 else:
  return HttpResponse("<p>Please log in to save recommendations.</p>")

####################################################
####################################################
####################################################
####################################################
############   REWRITTEN THUS FAR ##################
####################################################
####################################################
####################################################
####################################################

######################  Core Views  ####################################
import time

@login_required
def visuals(request):
 #Variable Setup
 u = request.user
 fatal_message = ""

 #TODO: Nora, you'll want to have some "loading page" that uses JavaScript
 #  to load the data via JSON. D3 makes this nifty easy. You'll want to make
 #  sure that a lab only can view THEIR OWN data (and eventually public data --
 #  but ignore this for now).

 return render(request, 'predictions_global.html', {
  "fatal_message": fatal_message,
 })

######################  Searching  #####################################
  #Rules:
  # 1.) A Lab can search only its data and public data.
  # 2.) Data that is returned can be sent to editing functions.

  #Verbose JSON Formats:
  # request.body ===
  #  query_list --> {[{u"field":u"FIELD", u"value":u"VALUE"}, ...]}

@login_required
def search(request, model="Data", params={}):
 from DRP.models import get_model_field_names

 u = request.user
 if request.method=="POST":
  try:

   #Get the appropriate data for the appropriate model.
   if model=="Data":
    #Pass the request on to data_transmit.
    return data_transmit(request)

   elif model=="Recommendation":
    #If applicable, only show reactions that are from seeds.
    seeded = params["seeded"]==True if "seeded" in params else False
    return recommendation_transmit(request, seeded=seeded)

  except:
   return HttpResponse("Woops! A problem occurred.")
 else:

  search_fields = get_model_field_names(both=True, unique_only=True, model=model)
  if model=="Data":
   #Collect the fields that will be displayed in the Search "Fields" tab.
   search_fields = get_model_field_names(both=True, unique_only=True)
   search_fields = [
   {"raw":"reactant", "verbose":"Reactant"},
   {"raw":"quantity", "verbose":"Quantity"},
   {"raw":"unit", "verbose":"Unit"},
   {"raw":"is_valid", "verbose":"Is Valid"},
   {"raw":"user", "verbose":"User"},
   {"raw":"public", "verbose":"Public"}] + search_fields
  elif model=="Recommendation":
   #Collect the fields that will be displayed in the Search "Fields" tab.
   search_fields = [
   {"raw":"reactant", "verbose":"Reactant"},
   {"raw":"quantity", "verbose":"Quantity"},
   {"raw":"unit", "verbose":"Unit"},
   {"raw":"assigned_user", "verbose":"Assigned User"}] + search_fields

   if "seeded" in params and params["seeded"]:
    model = "SeedRecommendation"
    search_fields += [{"raw":"seed", "verbose":"Seed is..."}]

  return render(request, 'search_global.html', {
   "search_fields": search_fields,
   "model": model,
  })



######################  CG Guide  ######################################
#Return a json object of the first ChemSpider result.
@login_required
@require_http_methods(["GET"])
def check_compound(request):
 u = request.user
 #Search through each available ChemSpider field.
 search_fields = [request.GET.get("CAS_ID"), request.GET.get("compound")]
 query = get_first_chemspider_entry(search_fields)
 if query:
  query_results = {
    "imageurl": query.imageurl,
    "commonName": query.commonname,
    "mv": query.molecularweight,
    "mf": query.mf,
   }
  response = json.dumps(query_results)
  return HttpResponse(response, content_type="application/json")
 return HttpResponse(1) #Return a code to signal the compound was not found.

#Ask for a confirmation whether a CG is correct or not.
@login_required
def compound_guide_form(request):
 u = request.user
 success = False
 lab_group = u.get_profile().lab_group
 if request.method == 'POST':
  #Bind the user's data and verify that it is legit.
  form = CompoundGuideForm(lab_group=lab_group, data=request.POST)
  #If all data is valid, save the entry.
  if form.is_valid():
   #Variable Setup.
   entry = form.save(commit=False) #Don't save to database, but make CompoundEntry.
   atoms = request.POST.get("atoms")
   mw = request.POST.get("mw")

   if not entry.custom:
    #Apply calculations to the compound.
    entry.save()
    update_compound_and_reactions(lab_group, entry)
    success = True #Used to display the ribbonMessage.
   elif atoms and mw.replace(".","",1).isdigit():
    entry.smiles = atoms
    entry.mw = mw
    entry.save()
    success = True #Used to display the ribbonMessage.
   else:
    return HttpResponse("4")
   #Clear the cached CG data.
   set_cache(lab_group, "COMPOUNDGUIDE", None)
   set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", None)

 else:
  #Submit a blank form if one was not just submitted.
  form = CompoundGuideForm()

 return render(request, 'compound_guide_form.html', {
   "form": form,
   "success": success,
 })

#Send/receive the compound guide:
@login_required
def compound_guide(request):
 u = request.user
 lab_group = u.get_profile().lab_group
 guide = list(get_lab_CG(lab_group))
 return render(request, 'compound_guide.html', {
  "guide": guide,
 })

@login_required
@require_http_methods(["POST"])
def compound_guide_entry(request):
 from DRP.models import get_lab_Data

 u = request.user
 lab_group = u.get_profile().lab_group
 entry_info = json.loads(request.body, "utf-8")
 query = get_lab_CG(lab_group).filter(compound=entry_info["compound"])
 if query.exists():
  return render(request, 'compound_guide_row.html', {
   "entry": query[0]
  })
 else:
  return HttpResponse("<p>No CG found. Please refresh page.</p>")

def change_Data_abbrev(lab_group, old_abbrev, new_abbrev):
 if old_abbrev=="":
  raise Exception("Abbrev cannot be an empty string.")

 lab_data = get_lab_Data(lab_group)
 for i in CONFIG.reactant_range():
  reactant = "reactant_{}".format(i)
  affected_data = lab_data.filter(Q((reactant, old_abbrev)))
  if affected_data.exists():
   affected_data.update(**{reactant:new_abbrev})


### EDIT ME
@login_required
@require_http_methods(["POST"])
def edit_CG_entry(request):
 from DRP.models import get_lab_Data

 u = request.user
 changesMade = request.POST

 #Get the Lab_Group's Compound Guide
 lab_group = u.get_profile().lab_group
 CG_data = get_lab_CG(lab_group)

 if changesMade["type"]=="del":
  try:
   lab_data = get_lab_Data(lab_group)
   #Delete each datum in the pid list.
   for pid in changesMade.getlist("pids[]"):
    #Mark any entry that uses this datum is now invalid.
    entry = CG_data.get(id=pid)
    affected_data = get_Data_with_abbrev(lab_data, entry.abbrev)
    affected_data.update(is_valid=False)

    #Now delete the CG entry.
    entry.delete()
  except Exception as e:
   print e
   return HttpResponse(1)
 elif changesMade["type"]=="edit":
  try:
   #Variable Setup
   field = changesMade["field"]
   new_val  = changesMade["newVal"]
   pid = changesMade["pid"]

   #Collect the datum to be changed and the old value.
   changed_entry = CG_data.get(id=pid)
   old_val = getattr(changed_entry, field)

   #Make sure the compound/abbrev isn't already being used.
   possible_entry = None
   if field=="compound":
    possible_entry = CG_data.filter(compound=new_val)
   elif field=="abbrev":
    possible_entry = CG_data.filter(abbrev=new_val)
   elif field=="CAS_ID" and new_val: #Empty CAS_IDs don't require queries.
    possible_entry = CG_data.filter(CAS_ID=new_val)
   if possible_entry!=None and possible_entry.exists():
    return HttpResponse("Already used!")

   #Make sure the new value does not invalidate the entry.
   dirty_data = model_to_dict(changed_entry)
   dirty_data[field] = new_val
   clean_data, errors = validate_CG(dirty_data, lab_group, editing_this=True)
   new_val = clean_data[field]
   if errors:
    raise Exception("Validation of datum failed.")

   #Commit the change to the Entry.
   setattr(changed_entry, field, new_val)
   changed_entry.save()
   #TODO:Remove this and make more "function" in style?
   set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", None)

   #Change the occurrences of the old abbrev to the new abbrev.
   if field=="abbrev":
    change_Data_abbrev(lab_group, old_val, new_val)

   #If no inorganic was found, ask the user for its identity.
   if changed_entry.compound_type=="Inorg":
    return HttpResponse("Inorganic not found!") #TODO: Add this.

   #Lookup fresh data from ChemSpider and RDKit
   update_compound(changed_entry)

  except Exception as e:
   print e
   return HttpResponse("Invalid!")

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
def data_form(request): #If no data is entered, stay on the current page.
 from DRP.models import get_lab_Data, get_model_field_names

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
   old_data_size = get_lab_Data_size(lab_group)
   clear_page_cache_of_index(lab_group, old_data_size)

   #Refresh the TOTALSIZE cache.
   set_cache(lab_group, "TOTALSIZE", old_data_size + 1)
   success = True #Used to display the ribbonMessage.
 else:
  try:
   lab_group = u.get_profile().lab_group
   model = request.GET["model"]
   pid = request.GET["pid"]

   #Either get initial fields from a Recommendation or a Data entry.
   if model=="rec":
    data = get_recommendations(lab_group).get(id=pid)
    init_fields = {field:getattr(data, field) for field in get_model_field_names(model="Recommendation")}
   else:
    data = get_lab_Data(lab_group).get(id=pid)
    init_fields = {field:getattr(data, field) for field in get_model_field_names()}

  except Exception as e:
   print e
   init_fields = {"leak":"No"}
  form = DataEntryForm(
   initial=init_fields
  )
 return render(request, 'data_form.html', {
  "form": form,
  "success": success,
 })

#Send/receive the data-entry form: #TODO: Merge with the field above.
@login_required
@require_http_methods(["GET"])
def transfer_rec(request):
 from DRP.models import get_model_field_names
 try:
  u = request.user
  lab_group = u.get_profile().lab_group
  pid = request.GET["pid"]
  rec = get_recommendations(lab_group).get(id=pid)

  initial_fields = {field:getattr(rec, field) for field in get_model_field_names(model="Recommendation")}
  form = DataEntryForm(
   initial=initial_fields
  )
  return render(request, 'data_form.html', {
   "form": form,
  })
 except:
  return HttpResponse("<p>Request failed!</p>")

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
@login_required
@require_http_methods(["POST"])
def upload_prompt(request):
 u = request.user
 model = request.POST.get("model")
 if model == "Data":
  return render(request, 'upload_form.html', {
   "model":model,
   "model_verbose":"Data",
  })
 elif model =="CompoundEntry":
  return render(request, 'upload_form.html', {
   "model":model,
   "model_verbose":"Compounds",
  })

 return HttpResponse("Request illegal!")

@login_required
def upload_CSV(request, model="Data"):
 u = request.user
 if request.method=="POST":
  #Variable Setup:
  lab_group = u.get_profile().lab_group
  fatal_message =""
  error_log=[]
  successes, fails, success_percent = 0, 0, 0

  #Attempt to get the data from the request
  print request.FILES #TODO DECLUTTER ME.
  return HttpResponse("Feature not added yet.")
  try:
   model=request.POST["model"]
   assert model in {"CompoundEntry", "Data"}
   uploaded_file = request.FILES["file"]
  except:
   return HttpResponse("5") #TODO: Add different types of responses.

  #More Variable Setup
  true_cols = get_model_field_names(model=model)
  file_fields = []
  auto_added_fields = {}
  list_field_nums = {field:1 for field in list_fields}
  current_row = 0
  current_col = 0

  return HttpResponse("Feature not added yet.")
  #Get the fields that are not required.
  not_required = get_not_required_fields(model=model) #TODO: WRITE ME.

  #Get the correct headings for the file.
  content = csv.reader(uploaded_file, delimiter=",")

  #Attempt to figure out the field order based on the headers in the file.
  for col in content[0]:
   try:
    field = get_related_field(col, model=model, for_upload=True)

    #Add the field number if it is missing.
    if field in list_fields:
     field += list_field_nums[field]
     list_field_nums[field] += 1

    true_cols.remove(field)
    file_fields.append(field) #Remember the order of fields in the file.

   except Exception as e:
    print e
    return HttpResponse("Valid headers were not found.")

  #Add the auto-add fields.
  auto_added_fields = set(true_fields)

  print field_fields
  print auto_added_fields

  #Make sure no vital fields are missed.
  for field in true_fields:
   assert field in not_required

#############################################################
  success_percent = 1 if not (fails+successes) else successes/(fails+successes)
  return render(request, 'upload_results.html', {
   "fatal_message": fatal_message, #Includes data and data_indexes.
   "error_log": error_log,
   "added_quantity": added_quantity,
   "error_quantity": error_quantity,
   "success_percent": success_percent,
  })
 return render(request, 'upload_form.html')

#TODO: Remo this. This is a legacy versino of the upload script for specifically data.
def upload_CSV_bak(request, model="Data"): ###Not re-read.
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
     "image_url",
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
  fatal_message = "<p>Please log in to upload data.</p>"
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


######################  Update Data ####################################
  #Rules:
  # 1.) A Lab can only delete data it owns.
  # 2.) Users can only modify their own Lab's data.

  #Verbose JSON Formats:
  # request.POST ===
  #  add_reactant --> {pid, group, reactant, quantity, unit}
  #  delete_reactant --> {pid, group}
  #  delete_Data --> {[PID_0, PID_1, ... , PID_N]}
  #  change_Data --> {PID, fieldChanged, newValue}

#Delete a reactant group from a datum.
def add_reactant(request):
 from DRP.models import get_lab_Data

 u = request.user
 if request.method=="POST" and u.is_authenticated():
  #Variable Setup
  reactantDict = {}
  lab_group = u.get_profile().lab_group
  lab_data = get_lab_Data(lab_group)

  #Gather the request and user info.
  try:
   group = request.POST["group"]
   pid = request.POST["pid"]
   reactantDict["reactant"] = request.POST["reactant"]
   reactantDict["quantity"] = request.POST["quantity"]
   reactantDict["unit"] = request.POST["unit"]
  except Exception as e:
   print e
   return HttpResponse(3)

  #Validate and the apply datum.
  try:
   #Do a basic validation of the reactant.
   datum = lab_data.get(id=pid)
   assert validate_name(reactantDict["reactant"], lab_group)
   assert reactantDict["unit"] in edit_choices["unitChoices"]
   assert quick_validation("quantity", reactantDict["quantity"])
  except Exception as e:
   return HttpResponse(4)

  #Actually add the fields to the datum.
  try:
   for entry in list_fields:
    setattr(datum, "{}_{}".format(entry, group), reactantDict[entry])
   datum.user = u
   datum.save()
  except:
   return HttpResponse(2)

  #Attempt to update/re-validate the full datum (but don't die on fail).
  try:
   update_reaction(datum, lab_group)
  except:
   pass
  return HttpResponse("0_close")
 else:
  try:
   group = request.GET["group"]
   pid = request.GET["pid"]
   return render(request, "add_reactant_form.html", {
    "group":group,
    "pid":pid,
    "unitChoices":edit_choices["unitChoices"],
   })
  except Exception as e:
   return HttpResponse("Illegal group specified.")

#Delete a reactant group from a datum.
@login_required
@require_http_methods(["POST"])
def delete_reactant(request):
 from DRP.models import get_lab_Data

 u = request.user
 lab_group = u.get_profile().lab_group
 lab_data = get_lab_Data(lab_group)
 try:
  group = int(request.POST["group"])
  pid = request.POST["pid"]

  #If the reactant is required, don't delete it.
  if group<=CONFIG.reactants_required:
   return HttpResponse("First two reactants required.")

  #Remove the reactant fields from the datum.
  datum = lab_data.get(id=pid)
  for field in list_fields:
   setattr(datum, "{}_{}".format(field, group), "")
  datum.user = u
  datum.save()

  #Attempt to update/re-validate the full datum (but don't die on fail).
  try:
   update_reaction(datum, lab_group)
  except:
   pass

  return HttpResponse(0)
 except Exception as e:
  return HttpResponse("Edit failed.")

@login_required
@require_http_methods(["POST"])
def delete_Data(request):
 from DRP.models import get_lab_Data

 u = request.user
 #Variable Setup
 lab_group = u.get_profile().lab_group
 deleteList = json.loads(request.body, "utf-8")
 lab_data = get_lab_Data(lab_group)

 #Find and delete data entries in a User's Lab.
 for pid in deleteList:
  try:
   lab_data.get(id=pid).delete()
  except:
   HttpResponse("<p>One or more selected data not found.</p>")

 #Finally, return a success code.
 return HttpResponse(0);

@login_required
@require_http_methods(["POST"])
def change_Recommendation(request):
 u = request.user
 try:
  #Variable Setup
  recommendations = get_recommendations(u.get_profile().lab_group)
  editLog = request.POST
  whitelist = {"notes"} #Only allow notes to be modified.

  #Read the editLog
  try:
   pid = editLog["pid"]
   fieldChanged = editLog["field"]
   newValue = editLog["newValue"]
   rec = recommendations.get(id=pid)

   #Verify that the fieldChanged is in the whitelist.
   assert fieldChanged in whitelist
  except:
   return HttpResponse("Datum not editable!")

  #Save the new value.
  setattr(rec, fieldChanged, newValue)
  rec.user = u
  rec.save()

  return HttpResponse(0)

 except:
  return HttpResponse("Edit unsuccessful...")

# Used to change fields in Data objects.
@login_required
@require_http_methods(["POST"])
def change_Data(request):
 from DRP.models import get_lab_Data

 #Fields that may be changed via this script.
 whitelist = set(get_model_field_names())

 u = request.user
 #Variable Setup
 lab_group = u.get_profile().lab_group
 editLog = request.POST
 lab_data = get_lab_Data(lab_group)

 #Get the Datum for the lab.
 try:
  pid = editLog["pid"]
  fieldChanged = editLog["field"]
  newValue = editLog["newValue"]
  datum = lab_data.get(id=pid)
  #Verify that the fieldChanged is in the whitelist.
  assert fieldChanged in whitelist
 except:
  return HttpResponse("Datum not editable!")

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
   return HttpResponse("{}".format(errors[errors.keys()[0]]))

  #Get the parsed value after cleaning.
  setattr(datum, fieldChanged, clean_data[fieldChanged])

  #Make the edit in the database.
  datum.user = u
  datum.is_valid = clean_data["is_valid"]
  datum.save()

  if fieldChanged=="ref":
   #Update the "ref" in any Data of which it is a duplicate.
   lab_data.filter(duplicate_of=oldValue).update(duplicate_of=newValue)
  elif fieldChanged[:8]=="reactant":
   #Update the atom information and any other information that may need updating.
   update_reaction(datum, lab_group)
  return HttpResponse(0)

 except Exception as e:
  print e
  return HttpResponse("Edit unsuccessful...")

