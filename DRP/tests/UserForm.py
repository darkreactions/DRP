# Set the Python path so that it has access to the Djano settings.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from django.contrib.auth.models import User

import django.forms as forms
from DRP.forms import *

from DRP.models import Data, Lab_Group, Lab_Member, CompoundEntry
from DRP.settings import ACCESS_CODE_MAX_LENGTH
from DRP.validation import validate_CG, full_validation
from DRP.validation import TYPE_CHOICES, UNIT_CHOICES
from DRP.validation import BOOL_CHOICES, PURITY_CHOICES, OUTCOME_CHOICES


#Testing if making a new user works 
def main():
  try:
    #Make a Recommendation
    new_user = UserForm(forms.ModelForm) 
    new_user.username = "new"
    new_user_Meta = new_user.Meta()
    new_user.save()

  except Exception as e:
    print "ERROR: {}".format(e)

if __name__=="__main__":
  main()



