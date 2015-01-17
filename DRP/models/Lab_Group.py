from django.db import models
from DRP.settings import ACCESS_CODE_MAX_LENGTH

class Lab_Group(models.Model):
  class Meta:
    app_label = "DRP"

  lab_title = models.CharField(max_length=200, unique=True, error_messages={'unique':"This name is already taken."})
  lab_address = models.CharField(max_length=200)
  lab_email = models.CharField(max_length=254) #Maximum length of email address


  def get_random_code(length = ACCESS_CODE_MAX_LENGTH):
    #Create a random alphanumeric code of specified length.
    import random, string
    options = string.letters + string.digits
    chars = [random.choice(options) for i in range(length)]
    return "".join(chars)

  access_code = models.CharField(max_length=ACCESS_CODE_MAX_LENGTH,
    default=get_random_code)

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

