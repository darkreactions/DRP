def migrateModelFields():
  def migrate(elem):
    from DRP.models import CompoundEntry
    from data_config import CONFIG

    for i in CONFIG.reactant_range():

      old_field = "reactant_{}".format(i)
      new_field = "reactant_fk_{}".format(i)

      abbrev = getattr(entry, old_field)

      comp = CompoundEntry.objects.filter(abbrev=abbrev).first()

      if not comp:
        comp = CompoundEntry()
        comp.abbrev = abbrev
        comp.compound = abbrev
        comp.lab_group = entry.lab_group
        comp.custom = True
        comp.save()

      setattr(entry, new_field, comp)

      entry.save()

  from DRP.models import Data

  database = Data.objects.all()

  for i, entry in enumerate(database):
    if (i%10==0):
      print "{}...".format(i)

    migrate(entry)

  print "FINISHED!"



def update_all_reactions(lab_group=None):
  from DRP.models import get_lab_Data, Data

  if lab_group:
    data = get_lab_Data(lab_group)
  else:
    data = Data.objects.all()

  for i, entry in enumerate(data):
    if (i%250==0): print "{}...".format(i+1)
    try:
      entry.refresh()
    except:
      print "--could not update reaction: {}".format(entry)

  print "Finished data validation."



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

  print "CG_calculations complete! ({} of {} entries changed)".format(success, cg.count())



def refresh_compound_guide(lab_group=None, verbose=False, debug=False, clear=False):
  #Either force-refresh all of the data or the data for a specific lab group.

  from DRP.models import *

  if lab_group:
    query = get_lab_CG(lab_group)
  else:
    query = CompoundEntry.objects.all()

  if clear:
    if debug: print "Clearing all CG_calculations..."
    CompoundEntry.objects.all().update(calculations=None)
    CompoundEntry.objects.all().update(calculations_failed=False)
    CG_calculations.objects.all().delete()

  #Actually perform the refresh/update.
  for i, compound in enumerate(query):
    try:
      if verbose and i%10==0: print "... {}.".format(i)
      update_compound(compound, debug=debug)
    except Exception as e:
      if debug: print "Could not update: {}\n\t".format(compound, e)


def recalculate_valid_data():
  from DRP.models import Data

  data = Data.objects.all()

  for i, datum in enumerate(data):
    if (i%250)==0: print "{}".format(i+1)
    try:
      datum.get_calculations_list()
      datum.is_valid = True
    except Exception as e:
      print e
      datum.is_valid = False
    print "{} {}".format(i, datum.is_valid)
    datum.save()

  print "{} Data valid.".format(Data.objects.filter(is_valid=True).count())
  print "{} Data invalid.".format(Data.objects.filter(is_valid=False).count())


def recalculate_usable_models():
  from DRP.models import ModelStats

  # Variable Setup
  all_models = ModelStats.objects.all()
  good = 0

  for model in all_models:
    usable = model.check_usability()
    if usable: good += 1

  print "{} of {} models usable!".format(good, all_models.count())
