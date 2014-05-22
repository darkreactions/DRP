import sys
from django.core.mail import send_mail
from DRP.retrievalFunctions import *
from DRP.errorReporting import print_error

######################  Email Functions  ########################### 
def email_admins(subject, message):
  try:
    send_mail("DRP: {}".format(subject), message, settings.EMAIL_HOST_USER, [settings.EMAIL_HOST_USER], fail_silently=False)
  except Exception as e:
    print_error("email_admins failed: {}\n".format(e))


def email_user(user, subject, message):
  try:
    send_mail("DRP: {}".format(subject), message, settings.EMAIL_HOST_USER, 
              [user.email], fail_silently=False)
  except:
    print_error("email_user failed: {}\n".format(user.id))


def email_lab(lab, subject, message, include_members=False):
  try:
    send_mail("DRP: {}".format(subject), message, settings.EMAIL_HOST_USER, 
              [lab.lab_email], fail_silently=False)
  except:
    print_error("email_lab failed: {} ({})\n".format(lab.id, lab.lab_title))

  #Email every lab member if the option is specified.
  if include_members:
    users = get_lab_users(lab_group)
    for user in users:
      email_user(user, "(Lab) {}".format(subject), message)


######################  Email Wrappers  ########################### 
def alert_about_new_lab(lab_group):
  email_body = "A new lab group has been registered:\n{}\n{}\n{}".format(lab_group.lab_title, lab_group.lab_address, lab_group.lab_email) 
  email_admins("New Lab Group Registered", email_body)
