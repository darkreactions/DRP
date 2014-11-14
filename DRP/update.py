
def update_data_calcs(only_invalid=True, debug=False, clear=False):
  """
  Re-calculate every data entry.
  """

  from models import Data, DataCalc
  from model_building.load_cg import get_cg
  from model_building.load_data import get_abbrev_map

  if clear:
    if debug: print "Clearing Data Calculations..."
    Data.objects.all().update(calculations=None)
    DataCalc.objects.all().delete()


  error, total = 0,0
  errors = []

  cg = get_cg()
  abbrev_map = get_abbrev_map()

  for i, datum in enumerate(Data.objects.all()):
    try:
      if debug and (i%25)==0: print "{}..".format(i)
      datum.get_calculations_dict(debug=debug, preloaded_cg=cg,
                                               preloaded_abbrev_map=abbrev_map)
    except Exception as e:
      if debug: print "Error Data {} (valid: {}): {}".format(i, datum.is_valid, e)
      errors.append(datum)
    finally:
      total += 1

  success = total-len(errors)
  print "Success: {} / {} ({}%)".format(success, total, success/float(total)*100)

