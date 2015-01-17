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



def update_all_reactions(lab_group):
  from DRP.models import get_lab_Data
  data = get_lab_Data(lab_group)
  for entry in data:
    try:
      entry.refresh()
    except:
      print "--could not update reaction: {}".format(entry)
  print "Finished data validation."



def perform_CG_calculations(only_missing=True, lab_group=None,
                            attempt_failed = True, verbose=False):

  from DRP.validation import clean_compound
  from DRP.models import create_CG_calcs_if_needed

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
        entry.compound = clean_compound(entry.compound)
        calc = create_CG_calcs_if_needed(entry.compound, entry.smiles, entry.compound_type)
        entry.calculations = calc
      except:
        entry.calculations_failed = True

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

