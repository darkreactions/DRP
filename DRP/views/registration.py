# # # # # # # # # # # # # # # # # # # 
 # # #  Registration Views   # # # #
# # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
from django.http import HttpResponse
from django.shortcuts import render

from DRP.retrievalFunctions import *
from DRP.emailFunctions import *

from DRP.forms import UserProfileForm, UserForm, LabForm
from django.contrib import auth

#Redirects user to the appropriate registration screen.
def registration_prompt(request):
 return render(request, "registration_cell.html", {})

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

    #Send a "confirmation" email to the new user.
    email_body = "This email confirms that you ({}) are now a registered member of the Dark Reaction Project under the following lab: {}".format(new_user.first_name, lab_group.lab_title)
    email_user(new_user, "User Registration Successful", email_body)

    #Politely log the user in!
    new_user = auth.authenticate(username = request.POST["username"],
     password = request.POST["password"])
    auth.login(request, new_user)
    #Mark that the user must now agree to the terms and conditions.
    return HttpResponse(1);
   else:
    return HttpResponse("Invalid Access Code!")
 else:
  form = [UserForm(), UserProfileForm()]
 return render(request, "user_registration_form.html", {
  "user_form": form[0],
  "profile_form": form[1],
 })

def lab_registration(request):
 if request.method=="POST":
  form = LabForm(data = request.POST)
  if form.is_valid():
   lab_group = form.save()
   access_code = lab_group.access_code
   
   #Send a "confirmation" email to the new lab email.
   email_body = "Thank you for joining the Dark Reaction Project!\n\nPlease continue by creating a \"user\" for your lab. Simply...\n\t1.) Keep this \"access code\" for your lab in a safe place: {}\n\t2.) Click \"Register\" and create a user using the access code above.\n\t3.) Start uploading data!\n\nWe wish you all the best,\nThe Dark Reaction Project Team".format(lab_group.access_code)

   email_lab(lab_group, "Lab Registration Successful", email_body)

   #Send the DRP Admins an email about the new Lab Registration.
   #TODO:Remove this when we scale to unmanageable quantities of labs.
   alert_about_new_lab(lab_group)
 
   return HttpResponse("<p>Registration Successful! Please check your email.</p>")
 else:
  form = LabForm()
 return render(request, "lab_registration_form.html", {
  "form": form,
 })
