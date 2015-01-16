# # # # # # # # # # # # # # # # # # #
 # # # # # # JSON Views  # # # # # #
# # # # # # # # # # # # # # # # # # #

#Necessary Imports:
from django.http import HttpResponse

from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
from django.shortcuts import render

from DRP.retrievalFunctions import *

import csv, json

@login_required
@require_http_methods(["POST"])
def download_CSV(request):
  from DRP.models import get_lab_Data, get_lab_CG, get_model_field_names

  #Variable Setup
  u = request.user
  lab_group = u.get_profile().lab_group

  #Get the info from the POST request.
  try:
    filters = request.POST.get("filters")
    if filters:
      filters = json.loads(filters)
    model = request.POST.get("model")
    assert model in {"Recommendation", "Data", "Saved", "CompoundEntry"}
  except Exception as e:
    return HttpResponse("Download request failed!")

  #File Setup
  file_name = "{}_{}".format(lab_group.lab_title, model).replace(" ", "_").lower()
  CSV_file = HttpResponse(content_type="text/csv")
  CSV_file["Content-Disposition"] = "attachment; filename={}.csv".format(file_name)
  result = csv.writer(CSV_file)

  #Write the headers to the file.
  if model=="Saved":
    model="Recommendation"
    saved_only = True
  else:
    saved_only = False

  verbose_headers = get_model_field_names(verbose=True, model=model)
  headers = get_model_field_names(verbose=False, model=model)

  #Modify the headers if needed and get/filter the data.
  if model=="Data":
    #Make sure the "Reference" is the first column.
    verbose_headers.remove("Reference")
    verbose_headers.insert(0, "Reference")
    headers.remove("ref")
    headers.insert(0, "ref")
    data = get_lab_Data(lab_group)

    if filters:
      data = filter_data(lab_group, filters)

  elif model=="Recommendation":
    print filters

    #If there are any filters, apply them.
    if filters:
      data = filter_recommendations(lab_group, filters)
    else:
      data = get_recommendations_by_date(lab_group)


  else: #if model=="CompoundEntry"
    verbose_headers.remove("Image URL")
    headers.remove("image_url")
    data = get_lab_CG(lab_group)

  try:
    #Write the data to the file.
    result.writerow(verbose_headers)
    for entry in data:
      result.writerow([getattr(entry, field).encode('ascii', errors='ignore') for field in headers])
    return CSV_file
  except Exception as e:
    print e
    return HttpResponse("Error preparing the file!")

@login_required
@require_http_methods(["POST"])
def download_prompt(request):
 u = request.user
 model = request.POST.get("model")
 print model
 if model in {"Data"}:
  return render(request, 'download_form.html', {
   "model":model,
   "model_verbose":"Data",
   "allow_filters":True
  })
 elif model in {"CompoundEntry","Saved","Recommendation"}:
  return render(request, 'download_form.html', {
   "model":model,
   "model_verbose":{
     "CompoundEntry":"Compounds",
     "Saved":"Saved",
     "Recommendation":"Recommendations"}[model],
   "allow_filters":True
  })

 return HttpResponse("Request illegal!")
