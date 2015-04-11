from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
from django.shortcuts import render

import json

# # # # # # # # # # # # # # # # # # #
  # # # # # # # # View Functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #

@login_required
@require_http_methods(["POST"])
def edit_recommendation(request, action):
 from DRP.models import get_recommendations

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

# # # # # # # # # # # # # # # # # # #
   # # # # # Sub-view Functions (eg, Javascript response views) # # # # # # # #
# # # # # # # # # # # # # # # # # # #
#Change the assigned_user of a recommendation.
@login_required
@require_http_methods(["POST"])
def assign_user_to_rec(request):
 from DRP.models import User, get_recommendations

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

######################  Searching  #####################################
  #Rules:
  # 1.) A Lab can search only its data and public data.
  # 2.) Data that is returned can be sent to editing functions.

  #Verbose JSON Formats:
  # request.body ===
  #  query_list --> {[{u"field":u"FIELD", u"value":u"VALUE"}, ...]}

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

######################  Database Functions  ############################
#Send/receive the data-entry form:
def data_form(request): #If no data is entered, stay on the current page.
 from DRP.models import get_lab_Data, get_model_field_names,get_recommendations
 from DRP.forms import DataEntryForm
 from DRP.retrievalFunctions import get_lab_Data_size
 from DRP.cacheFunctions import set_cache

 u = request.user
 success = False
 if request.method == 'POST' and u.is_authenticated():
  #Bind the user's data and verify that it is legit.
  form = DataEntryForm(user=u, data=request.POST)
  if form.is_valid():
   print "VALID"
   #If all data is valid, save the entry.
   form.save()
   lab_group = u.get_profile().lab_group

   #Clear the cache of the last page.
   old_data_size = get_lab_Data_size(lab_group)

   #Refresh the TOTALSIZE cache.
   set_cache(lab_group, "TOTALSIZE", old_data_size + 1)
   success = True #Used to display the ribbonMessage.
 else:
  print "INVALID"
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

  except:
   init_fields = {"leak":"No"}

  form = DataEntryForm(
   initial=init_fields
  )

 return render(request, 'data_form.html', {
  "form": form,
  "success": success,
 })

 ##################  Helper Functions ###############################


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
 from DRP.models import get_recommendations
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

