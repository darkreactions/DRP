from django.core.mail import send_mail
from DRP.retrievalFunctions import *

######################  Email Functions  ###########################
def alert_about_new_lab(lab_group):
  email_body = "A new Lab Group has been registered:"
  email_body += "\n{}\n{}\n{}".format(lab_group.lab_title, lab_group.lab_address, lab_group.lab_email) 
  send_mail("DRP: New Lab Group Registered", email_body, settings.EMAIL_HOST_USER, [settings.EMAIL_HOST_USER], fail_silently=False)

def email_user(user, subject, message):
  send_mail("DRP: {}".format(subject), message, settings.EMAIL_HOST_USER, [user.email], fail_silently=False)

def email_lab(lab, subject, message, only_head=True):
  send_mail("DRP: {}".format(subject), message, settings.EMAIL_HOST_USER, [lab_group.lab_email], fail_silently=False)

  #Email every lab member if the option is specified.
  if not only_head:
    users = get_lab_users(lab_group)
    for user in users:
      email_user(user, "(Lab) {}".format(subject), message)

