#!/usr/local/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # #  Compound Calculate 'n Store  Worker Process  # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


import django.db
from DRP.models import CompoundEntry, clean_compound, create_CG_calcs_if_needed

from DRP.emailFunctions import email_admins
from DRP.logPrinting import print_error, print_log

#An independent worker process for calculating compound properties.
def compound_calc_worker(compound_id, debug=False):
  print_log("Compound Calc: {}".format(compound_id))
  try:
    #Restart the database connection for this new process.
    django.db.close_connection()

    #Get the object from the database (assuming already-valid compounds).
    entry = CompoundEntry.objects.get(id=compound_id)

    try:
      entry.compound = clean_compound(entry.compound)
      calc = create_CG_calcs_if_needed(entry.compound, entry.smiles, entry.compound_type)
      entry.calculations = calc
      entry.calculations_failed = False
      status = "Passed"
    except Exception as e:
      status = "Failed"
      entry.calculations_failed = True

    #Save the calculation and mark it in the log.
    entry.save()
    if debug: print "#{} -- {}".format(compound_id, status)
  except Exception as e:
    #Log any errors that might have occurred.
    print_error("Compound ID Failed ({})".format(compound_id), details=e)
    email_admins("Fatal Failure: Compound Calc Worker",
                 "Compound ID \"{}\" failed auto-calculations.".format(compound_id))

if __name__ == "__main__":
  if not (1 < len(sys.argv) < 4):
    print "You probably want to let the site handle this..."
    print "USAGE: python ./calculate_compound_properties.py compound_id [debug]"
  else:
    debug = "debug" in sys.argv
    compound_calc_worker(sys.argv[1], debug=debug)



