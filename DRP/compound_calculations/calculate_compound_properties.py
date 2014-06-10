#!/usr/local/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # #  Compound Calculate 'n Store  Worker Process  # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
import sys, os

#Grab the Django settings if they aren't already set.
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

import django.db
from DRP.settings import BASE_DIR, LOG_DIR
from DRP.models import CompoundEntry, clean_compound, create_CG_calcs_if_needed

from DRP.emailFunctions import email_admins
from DRP.logPrinting import print_error, print_log

#An independent worker process for calculating compound properties.
def compound_calc_worker(compound_id):
  print_log("Compound Calc: {}".format(compound_id))
  try:
    #Restart the database connection for this new process.
    django.db.close_connection()

    try:
      #Get the object from the database (assuming already-valid compounds).
      entry = CompoundEntry.objects.get(id=compound_id)
      entry.compound = clean_compound(entry.compound)
      calc = create_CG_calcs_if_needed(entry.compound, entry.smiles, entry.compound_type)
      entry.calculations = calc
      status = "Passed"
    except:
      status = "Failed"
      entry.calculations_failed = True

    #Save the calculation and mark it in the log.
    entry.save
    print_log("Compound Calc: {} ({})".format(compound_id, status))
  except Exception as e:
    #Log any errors that might have occurred.
    print_error("Compound ID: {}".format(compound_id), details=e)
    email_admins("Fatal Failure: Compound Calc Worker",
                 "Compound ID \"{}\" failed auto-calculations.".format(compound_id))

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print "You probably want to let the site handle this..."
    print "USAGE: python ./this_script.py compound_id"
  else:
    compound_calc_worker(sys.argv[1])
    


