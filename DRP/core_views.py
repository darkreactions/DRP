from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
from django.shortcuts import render

import json
import string

from data_config import CONFIG

# # # # # # # # # # # # # # # # # # #
  # # # # # # # # Data and Page Helper Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # #
  # # # # # # # # View Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #

@login_required
def recommend(request, page_request=None):
  from DRP.pagifier import get_page_link_format
  from DRP.retrievalFunctions import get_recommendations_by_date

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
    total_pages = 1

  return render(request, 'global_page.html', {
  "template":"recommendations",
  "recommendations": recommendations,
  "fatal_message": fatal_message,
  "total_data_size":total_data_size,
  "counter_offset":recs_per_page*(page-1),
  "page_info": {
                    "data_per_page":recs_per_page,
                    "current_page":page,
                    "page_links":get_page_link_format(page, total_pages),
                    }
 })


@login_required
@require_http_methods(["POST"])
def recommendation_transmit(request, seeded=False):
 from DRP.pagifier import get_page_link_format
 from DRP.retrievalFunctions import get_recommendations_by_date, filter_recommendations

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
   "page_info": {
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
 from DRP.retrievalFunctions import get_recommendations

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
 from DRP.retrievalFunctions import get_recommendations
 from DRP.models import User

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
  import random
  def get_random_unranked_reaction_or_none(seed=None):
    def get_unranked_reactions(seed=None):
      from DRP.models import RankedReactionList
      unranked = RankedReactionList.objects.filter(ranker=None)
      #If a seed is specified, apply it to the filter.
      if seed:
        #Convert any lists to strings to allow filtering.
        seed = json.dumps(seed) if type(seed)==list else seed
        unranked = unranked.filter(seed=seed)
      return unranked

    unranked_rxns = get_unranked_reactions(seed=seed)
    if unranked_rxns.exists():
      random_index = random.randrange(unranked_rxns.count())
      return unranked_rxns[random_index]
    return None

  #Variable Setup:
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
 from DRP.retrievalFunctions import get_recommendations
 from DRP.models import User

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

@login_required
def visuals(request):
 #Variable Setup
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


  search_sub_fields = [
    {"raw":"abbrev", "verbose":"Abbrev"},
    {"raw":"compound", "verbose":"Compound"},
    {"raw":"compound_type", "verbose":"Type"}
  ]

  return render(request, 'search_global.html', {
   "search_fields": search_fields,
   "search_sub_fields": search_sub_fields,
   "model": model,
  })



######################  CG Guide  ######################################
#Return a json object of the first ChemSpider result.
@login_required
@require_http_methods(["GET"])
def check_compound(request):
  from DRP.chemspider import search_chemspider

  #Search through each available ChemSpider field.
  search_fields = [request.GET.get("CAS_ID"), request.GET.get("compound")]
  query = search_chemspider(search_fields)

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
 from DRP.forms import CompoundGuideForm
 from DRP.cacheFunctions import set_cache

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


@login_required
@require_http_methods(["POST"])
def compound_guide_entry(request):
 from DRP.models import get_lab_CG

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


#TODO: HOPEFULLY WILL MOSTLY DEPRECATE ON DB MIGRATION.
@login_required
@require_http_methods(["POST"])
def edit_CG_entry(request):
 from DRP.models import get_lab_CG, get_lab_Data, get_Data_with_Compound
 from DRP.cacheFunctions import set_cache
 from DRP.validation import validate_CG

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
    affected_data = get_Data_with_Compound(lab_data, entry.abbrev)
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
   clean_data, errors = validate_CG(changed_entry, lab_group, editing_this=True)
   new_val = clean_data[field]
   if errors:
    raise Exception("Validation of datum failed.")

   #Commit the change to the Entry.
   setattr(changed_entry, field, new_val)
   changed_entry.save()
   #TODO:Remove this and make more "function" in style?
   set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", None)

   #If no inorganic was found, ask the user for its identity.
   if changed_entry.compound_type=="Inorg":
    return HttpResponse("Inorganic not found!") #TODO: Add this.

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
 from DRP.forms import DataEntryForm
 from DRP.retrievalFunctions import get_lab_Data_size, get_recommendations
 from DRP.cacheFunctions import set_cache

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

    # Replace compound entries with an `abbrev`.
    for field, val in init_fields.items():
      if "reactant" in field:
        init_fields[field] = getattr(val, "abbrev")

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
 from DRP.retrievalFunctions import get_recommendations
 from DRP.forms import DataEntryForm
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
@require_http_methods(["POST"])
def upload_CSV(request, model="Data"):
  return HttpResponse("Feature not added yet.")

######################  Update Data ####################################

@login_required
@require_http_methods(["POST"])
def change_Recommendation(request):
 from DRP.retrievalFunctions import get_recommendations
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

