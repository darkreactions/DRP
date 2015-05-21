# # # # # # # # # # # # # # # # # # #
 # # # # # # JSON Views  # # # # # #
# # # # # # # # # # # # # # # # # # #

#Necessary Imports:

from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods


######################  Data Transmit ##################################
#Send the CG name pairs to the client.
@login_required
@require_http_methods(["GET"])
def send_CG_names(request):
  from django.http import HttpResponse
  from DRP.models import collect_CG_name_pairs
  import json

  u = request.user
  lab_group = u.get_profile().lab_group
  name_pairs = collect_CG_name_pairs(lab_group)
  return HttpResponse(json.dumps(name_pairs), mimetype="application/json")

