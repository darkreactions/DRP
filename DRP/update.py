
def update_data_calcs(only_invalid=True):
  """
  Re-calculate every data entry.
  """

  from DRP.models import Data


  error, total = 0,0
  errors = []

  for datum in Data.objects.all():
    try:
      datum.get_calculations_dict(force_recalculate=True)
    except:
      errors.append(datum)
    finally:
      total += 1

  success = total-len(errors)
  print "Complete: {} / {} ({}%)".format(success, total, success/float(total)*100)


def update_compound_guide():


