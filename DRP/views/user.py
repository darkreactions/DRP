# # # # # # # # # # # # # # # # # # # 
 # # # User Views and Functions  # #
# # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from django.contrib import auth
from django.db.models import Q

from django.contrib.auth.decorators import login_required

from DRP.retrievalFunctions import *
from DRP.emailFunctions import *

from DRP.forms import UserForm
from DRP.views.license_agreement import user_license_is_valid

######################  User Auth ######################################
def change_password(request):
 error=False
 if request.method == "POST":
  try:
   email = request.POST.get("email")
   username = request.POST.get("username")
   last_name = request.POST.get("lastName")
   user = User.objects.filter(Q(email=email)|Q(username=username), Q(last_name=last_name))[0]
   #Change the user's password and send them an email.
   randomize_password(user)
   return HttpResponse("A new password has been emailed to you.")
  except:
   #If no user is found given the credentials, tell the user.
   error=True
 return render(request, "change_password_form.html", {
  "error":error
 })


#Given a user, change their password and email them the new password.
def randomize_password(user):
 new_pass = get_random_code(15) #Generate a random password for the user.
 user.password = make_password(new_pass) #Hash the password. 
 user.save()
 email_body = "Hello {},\n\n According to our records, you just requested a password change. We have changed your account information as follows:\nUsername: {}\nPassword: {}".format(user.first_name, user.username, new_pass)
 email_user(user, "Password Change Request", email_body)


def user_login(request):
 login_fail = False #The user hasn't logged in yet...

 if request.method == "POST":
  username = request.POST.get("username", "")
  password = request.POST.get("password", "")
  redirect_URL = request.POST.get("next", "")
  user = auth.authenticate(username=username, password=password)

  if user is not None and user.is_active:
   auth.login(request, user)

   if not user_license_is_valid(user):
    return HttpResponse(1);

   else:
     return HttpResponse("Logged in successfully! <div class=reloadImmediately redirect=\""+redirect_URL+"\"></div>");
  else:
   login_fail = True #The login info is not correct.
 return render(request, "login_form.html", {
  "login_fail": login_fail,
 })

@login_required
def user_logout(request):
 auth.logout(request)
 return HttpResponse("OK")

@login_required
def user_update(request):
 u = request.user
 if request.method == "POST":
  form = UserForm(request.POST, instance=u)
  if form.is_valid():
   form.save()
   return HttpResponse("Update Successful!")
 else:
  form = UserForm(instance=u)
 return render(request, "user_update_form.html", {
  "form": form,
 })
