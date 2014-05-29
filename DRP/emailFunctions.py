import sys
from django.core.mail import send_mail
from DRP.retrievalFunctions import *
from DRP.logPrinting import print_error
from DRP.settings import DEBUG

"""
The Email Functions below are designed to make emailing the admins (us),
the lab heads, and the lab members easier. Each takes a subject and a 
message (and a lab/user if applicable), and sends them an email. If a
critical error occurs (which probably means the darkreactions email/pass
is incorrect), an error is printed. Each automatic "worker" process that
is added should include emails to admins to indicate failures.

Note that in DEBUG-mode (indicated in settings.py) email errors will not
be logged -- since emailing should never work if you have the password
correctly set to "SecurePassword" as you should while developing.
"""

######################  Email Functions  ########################### 
def email_admins(subject, message, include_managers=False):
  try:
    if include_managers:
      recipients = list(set(settings.ADMIN_EMAILS + settings.MANAGER_EMAILS))
    else:
      recipients = settings.ADMIN_EMAILS

    send_mail("DRP: {}".format(subject), message, settings.EMAIL_HOST_USER, recipients, fail_silently=False)
  except Exception as e:
    if not DEBUG:
      print_error("email_admins failed: {}\n".format(e))


def email_user(user, subject, message):
  try:
    send_mail("DRP: {}".format(subject), message, settings.EMAIL_HOST_USER, 
              [user.email], fail_silently=False)
  except:
    if not DEBUG:
      print_error("email_user failed: {}\n".format(user.id))


def email_lab(lab, subject, message, include_members=False):
  try:
    send_mail("DRP: {}".format(subject), message, settings.EMAIL_HOST_USER, 
              [lab.lab_email], fail_silently=False)
  except:
    if not DEBUG:
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
