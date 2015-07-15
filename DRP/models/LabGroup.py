from django.db import models
from django.contrib import auth
import random, string

# Helper function for creating random access codes for Lab Groups.
def get_random_code():
  #Create a random alphanumeric code of specified length.
  options = string.letters + string.digits
  chars = [random.choice(options) for i in xrange(ACCESS_CODE_LENGTH)]
  return "".join(chars)

class LabGroup(models.Model):
  class Meta:
    app_label = "DRP"

  title = models.CharField(max_length=200, unique=True, error_messages={'unique':"This name is already taken."})
  address = models.CharField(max_length=200)
  email = models.CharField(max_length=254) #Maximum length of email address
                                 default='')
  access_code = models.CharField(max_length=128)
  legacy_access_code = models.CharField(max_length=20)
  users = models.ManyToManyField(user)

  def __unicode__(self):
    return self.lab_title

def get_Lab_Group(query):
  try:
    if type(query)==Lab_Group:
      return query
    else:
      return Lab_Group.objects.filter(lab_title=query).get()
  except:
    message = "Could not find Lab_Group with lab_title: {}".format(query)
    raise Exception(message)

