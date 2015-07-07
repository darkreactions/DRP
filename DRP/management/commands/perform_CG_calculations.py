from django.core.management.base import BaseCommand
from optparse import make_option

def perform_CG_calculations(only_missing=True, lab_group=None,
                            attempt_failed = True, verbose=False):

  from DRP.models import CompoundEntry

  #Variable Setup
  success = 0
  i = 0

  cg = CompoundEntry.objects.all()
  if only_missing:
    cg = cg.filter(calculations=None)
  if lab_group:
    cg = cg.filter(lab_group=lab_group)
  if not attempt_failed:
    cg = cg.filter(calculations_failed=False)

  for entry in cg:
    try:
      if verbose:
        i+=1
        if i%5==0: print "... {}.".format(i)

      try:
        entry.create_CG_calcs_if_needed()

      except Exception as e:
        entry.calculations_failed = True
        print e

      entry.save
      success += 1
    except Exception as e:
      print "CG_calculation construction failed: {}".format(entry.compound)
      print "ERROR: {}".format(e)

class Command(BaseCommand):
  option_list = BaseCommand.option_list + (
  make_option('--only-missing',
       action='store_true',
       dest='only_missing',
       default=False,
       help='Perform calculations on only the missing data.'),
  ) + (
  make_option('--attempt-failed',
       action='store_true',
       dest='attempt_failed',
       default=True,
       help='Perform calculations on previously failed calculations.'),
  )

  def handle(self, *args, **options):
   self.stdout.write("Calculations started!")

   perform_CG_calculations(only_missing=options["only_missing"],
                           attempt_failed=options["attempt_failed"],
                           verbose=True)

   self.stdout.write("CG_calculations construction complete!")

