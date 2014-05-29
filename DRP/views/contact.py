# # # # # # # # # # # # # # # # # # # 
 # # # #  Contact Views  # # # # # 
# # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
from django.http import HttpResponse
from django.views.decorators.http import require_http_methods
import datetime

from DRP.retrievalFunctions import *
from DRP.emailFunctions import *

@require_http_methods(["POST"])
def contact_form(request):
  try:
    #Get the contact form information.
    email = request.POST["email"]
    content = request.POST["content"]
    subject = request.POST["subject"]

    #Make sure every field is filled, lest a fail should occur.
    if not email or not content or not subject:
      return HttpResponse("<p>Please fill out the entirety of the form!</p>")

    #Get the current datetime to help with comparison with log files.
    time = datetime.datetime.now()

    #Get the user's information if possible.
    if request.user.is_authenticated():
      u = request.user
      user_info = "{} (id: {})".format(u.id, u.username)
    else:
      user_info = "Anonymous User"

    body = "SUBJECT:{}\nEMAIL:{}\nUSER:{}\n\nCONTENT:{}\n\n({})".format(subject,
                                                                email, 
                                                                user_info, 
                                                                content,
                                                                time) 

    email_admins("Contact Form ({})".format(subject), body, include_managers=True) 

    return HttpResponse("<p>Thank you for contributing to the Dark Reaction Project. We'll get back to you soon.</p>")

  except Exception as e:
    print e
    return HttpResponse("<p>Oops... Something went wrong. Sorry about that!</p>")
 

