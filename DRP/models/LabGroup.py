'''A module containing only the LabGroup class'''
from django.db import models
from django.contrib import auth
from django.contrib.auth.models import User

class LabGroup(models.Model):
  '''A class for describing a collection of scientists belonging to the same group.'''
  class Meta:
    app_label = "DRP"

  title = models.CharField(max_length=200, unique=True, error_messages={'unique':"This name is already taken."})
  address = models.CharField(max_length=200)
  email = models.CharField(max_length=254, #Maximum length of email address
                                 default='')
  access_code = models.CharField(max_length=128)
  '''An access code to allow members to \'prove\' that they are permitted to join the lab. Normally held by a laboratory administrator.'''
  legacy_access_code = models.CharField(max_length=20)
  '''An older version of the access code. Made a part of this model for legacy support'''
  users = models.ManyToManyField(User, blank=True)

  def __unicode__(self):
    return self.title

def get_Lab_Group(query):
  try:
    if type(query)==Lab_Group:
      return query
    else:
      return Lab_Group.objects.filter(lab_title=query).get()
  except:
    message = "Could not find Lab_Group with lab_title: {}".format(query)
    raise Exception(message)

