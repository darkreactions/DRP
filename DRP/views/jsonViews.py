# # # # # # # # # # # # # # # # # # # 
 # # # # # # JSON Views  # # # # # #
# # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
from django.http import HttpResponse

from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods

from DRP.retrievalFunctions import *
from DRP.models import convert_QuerySet_to_list
######################  Data Transmit ##################################
#Send the CG name pairs to the client.
@login_required
@require_http_methods(["GET"])
def send_CG_names(request):
 u = request.user
 lab_group = u.get_profile().lab_group
 name_pairs = collect_CG_name_pairs(lab_group, overwrite=False)
 return HttpResponse(json.dumps(name_pairs), mimetype="application/json")

#TODO: Nora, this is where I'd throw the functions that you want to _send_ 
#      data. This way, we can keep everything organized.

#Just because I'm not sure if Brian's code will interact well with the database, I am going to 
#steal from DRP/download.py 's code to turn data into a CSV, and then pass to Brian's code for the vis
@login_required
def pull_CSV(lab_group)
  data = get_lab_Data(lab_group) 
  data_as_list = convert_QuerySet_to_list(data) 
  expanded_csv = [data_as_list[i] + data[i].calculations.make_json() for i in xrange(len(data_as_list))] 
  return expanded_csv #simply returns a list 	

# Get info from POST request (**must ask Casey what info this refers to (some options about the data?)
# Not sure how the filter works either, nor what the "model" is doing in the "assert" statement  
 
	

#File set_up (trying to adapt DRP/download.py code to create a CSV file, but not make it a download)
# Function to parse DataCalc's JSON text as CSV


#Function to parse CSV data into correct format for vis (place in explore view?) 



