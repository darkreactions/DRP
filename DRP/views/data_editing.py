from django.http import HttpResponse
from django.shortcuts import render

from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods


from DRP.validation import edit_choices

  #Rules:
  # 1.) A Lab can only delete data it owns.
  # 2.) Users can only modify their own Lab's data.

  #Verbose JSON Formats:
  # request.POST ===
  #  add_reactant --> {pid, group, reactant, quantity, unit}
  #  delete_reactant --> {pid, group}
  #  delete_Data --> {[PID_0, PID_1, ... , PID_N]}
  #  change_Data --> {PID, fieldChanged, newValue}


#TODO: Rewrite this.

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
      datum.refresh()
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
      print e
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
   datum.refresh()
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
   datum.refresh()
  return HttpResponse(0)

 except Exception as e:
  print e
  return HttpResponse("Edit unsuccessful...")

