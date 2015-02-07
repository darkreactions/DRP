# Set the Python path so that it has access to the Django settings.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def main():
  def _cleaner(data):
    headers = data.pop(0)
    matrix = [datum.get_calculations_list() for datum in data]

    mapper = {
      "yes":1,
      "no":0,
      "?":-1,
    }

    matrix = [[mapper[e] if e in mapper else e for e in row] for row in matrix]

    return [headers] + matrix


  def _stripper(splits, headers):
    cut = 20
    headers = headers[cut:]
    new_splits = {key:[row[cut:] for row in data] for key, data in splits.items()}
    return (new_splits, headers)

  def sklearn_test(data):
    model = ModelStats()
    model.construct("Pipeline Test",
                    data,
                    force = True,
                    tool="random forest",
                    library="sklearn",
                    preprocessor=_cleaner,
                    postprocessor=_stripper,
                    debug=True)
    model.summary()

  def weka_test(data):
    model = ModelStats()
    model.construct("Pipeline Test",
                    data,
                    force = True,
                    preprocessor=_cleaner,
                    postprocessor=_stripper,
                    library="weka",
                    tool="svc",
                    debug=True)
    model.summary()

  from DRP.models import ModelStats
  from DRP.retrievalFunctions import get_valid_data
  from DRP.model_building.rxn_calculator import headers

  tests = [sklearn_test, weka_test]
  errors = 0

  subset_size = 50
  data = [headers] + list(get_valid_data()[:subset_size])

  for function in tests:
    try:
      function(data)
    except:
      errors += 1

  print "______"*10
  print "{}/{} successes".format(len(tests)-errors, len(tests))


if __name__=="__main__":
  main()
