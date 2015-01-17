def get_Data_with_abbrev(lab_data, abbrev):
  if type(abbrev)==CompoundEntry:
    abbrev = abbrev.abbrev #Meta...

  Q_list = [Q(("reactant_{}".format(i),abbrev)) for i in CONFIG.reactant_range()]
  return lab_data.filter(reduce(operator.or_, Q_list))


def update_reactions_with_compound(lab_group, compound):
  #TODO: DEPRECATE ME
  #Update the individual "atom" records on each reaction.

  from DRP.models import get_Data_with_abbrev

  lab_data = get_lab_Data(lab_group)
  changed_reactions = get_Data_with_abbrev(lab_data, compound)
  for reaction in changed_reactions:
    try:
      reaction.refresh()
    except Exception as e:
      print "update_reactions_with_compound failed: {}".format(e)


#Update the compound by reloading the ChemSpider search data.
def update_compound(entry, debug=False):
  from DRP.fileFunctions import createDirIfNecessary
  from DRP.chemspider import search_chemspider
  from subprocess import Popen
  from settings import LOG_DIR, BASE_DIR, MODEL_DIR

  try:
    if not entry.custom: #Only update compounds that are not custom.
      #Get the most up-to-date ChemSpider info for a given CAS/compound.
      query = search_chemspider([entry.CAS_ID, entry.compound])
      if query:
        #Update the entry.
        entry.image_url, entry.smiles, entry.mw = query.imageurl, query.smiles, query.molecularweight
        perform_calcs = True

      else:
        perform_calcs = False
        if debug:
          print "Found legacy entry (should be `custom`): {}".format(entry.compound)

    else:
        perform_calcs = False
        entry.calculations = None
        entry.image_url, entry.smiles, entry.mw = "","",""
    entry.save()

    if perform_calcs:

      #Start a new compound-calc worker to determine compound properties.
      comp_log_dir = LOG_DIR+"/compound_calculations"
      createDirIfNecessary(comp_log_dir)

      err_log = open(comp_log_dir+"/error.log","a")
      act_log = open(comp_log_dir+"/process.log","a")
      worker_script = BASE_DIR+"/DRP/compound_calculations/calculate_compound_properties.py"
      command = "python {} {}".format(worker_script, entry.id)
      #Log to the files above and make the worker independent of the parent process.
      Popen(command.split(), stdout=act_log, stderr=err_log, close_fds=True)

  except Exception as e:
    entry.calculations_failed = False
    entry.save()
    raise Exception("Compound update ({}) failed: {}".format(entry, e))


