"""
These Email below are designed to make emailing the admins (us),
the lab heads, and the lab members easier. Each takes a subject and a 
message (and a lab/user if applicable), and sends them an email. If a
critical error occurs (which probably means the darkreactions email/pass
is incorrect), an uncaptured exception occurs. Each automatic "worker" process that
is added should include emails to admins to indicate failures.
"""

import sys
from django.core.mail import send_mail
from DRP import settings


class Email(object):
  """The base email class, sends email to a specified recipient from a specified sender""" 

  def __init__(self, subject, messageFrame, messageVars={}, to=[], sender=settings.EMAIL_HOST_USER):
    self.subject = subject
    self.message = messageFrame.format(messageVars)
    self.recipients = to
    self.sender = sender

  def send(self):
    return send_mail(subject, message, self.sender, self.recipients, fail_silently=FALSE)
  
class EmailToAdmins(Email):
  """Sends email specifically to administrators, with a specific flag for managers"""

  def __init__(self, subject, messageFrame, messageVars={}, includeManagers=False, sender=settings.EMAIL_HOST_USER):
    if includeManagers:
      super(EmailToAdmins, self).__init__(subject, messageFrame, messageVars, sender, settings.ADMIN_EMAILS)
    else:
      super(EmailToAdmins, self).__init__(subject, messageFrame, messageVars, sender, settings.ADMIN_EMAILS + settings.MANAGER_EMAILS)

#class EmailToLab(Email):
#
#  def __init__(self, subject, messageFrame, messageVars={}, includeMembers=False):
#    if includeMembers:
#      
#      super(EmailToLab, self).__init__(subject, messageFrame, messageVars)


######################  Email Wrappers  ########################### 
#def alert_about_new_lab(lab_group):
#  email_body = "A new lab group has been registered:\n{}\n{}\n{}".format(lab_group.lab_title, lab_group.lab_address, lab_group.lab_email) 
#  email_admins("New Lab Group Registered", email_body)
