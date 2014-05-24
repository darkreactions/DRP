# # # # # # # # # # # # # # # # # # # 
 # # # License Views and Functions  # #
# # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
from django.http import HttpResponse
from django.shortcuts import render

from django.contrib.auth.decorators import login_required

from DRP.retrievalFunctions import *
from DRP.emailFunctions import *

from DRP.data_config import CONFIG
import datetime

#Return whether the user_license is valid (True) or invalid/missing (False)
def user_license_is_valid(user):
 try:
  return user.get_profile().license_agreement_date_dt > CONFIG.current_license_date
 except:
  #Assume that if the query fails, the user is not licensed.
  return False

@login_required
def get_user_license_agreement(request):
 u = request.user
 #Indicate whether the user needs to agree to updated terms or sign the terms initially.
 if u.get_profile().license_agreement_date_dt:
  if not user_license_is_valid(u):
   license_changed = True
 else:
  license_changed = False
  
 return render(request, 'user_license_form.html', {
  "license_changed": license_changed,
  "license_file": CONFIG.current_license_file,
  "license_date": CONFIG.current_license_date, 
 })

def update_user_license_agreement(request):
 u = request.user 
 if request.method=="POST":
  try:
   u.get_profile().update_license()
   email_body = "This is a receipt for your records to indicate that you accepted our new License Agreement. If you are receiving this in error, please contact us immediately."
   email_user(user, "License Agreement Confirmed", email_body)
   return HttpResponse(
    "<p>You're all up-to-date!</p>" +
    "<div class=\"button refreshButton\">Explore</div>"
   )
  except:
   return HttpResponse("<p>Your request could not be completed. Please try again.</p>")
 else:
  return HttpResponse("<p>Please click the \"I Agree\" to accept the Terms and Conditions.</p>")
