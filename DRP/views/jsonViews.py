# # # # # # # # # # # # # # # # # # # 
 # # # # # # JSON Views  # # # # # #
# # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
from django.http import HttpResponse

from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods

from DRP.retrievalFunctions import *

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


