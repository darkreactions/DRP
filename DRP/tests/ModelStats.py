# Set the Python path so that it has access to the Django settings.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def main():
  def _cleaner(data):
    cut = 20

    headers = data.pop(0)[cut:]
    matrix = [datum.get_calculations_list()[cut:] for datum in data]

    mapper = {
      "yes":1,
      "no":0,
      "?":-1,
    }

    matrix = [[mapper[e] if e in mapper else e for e in row] for row in matrix]

    return [headers] + matrix

  from DRP.models import ModelStats
  from DRP.retrievalFunctions import get_valid_data
  from DRP.model_building.rxn_calculator import headers

  subset_size = 200
  data = [headers] + list(get_valid_data()[:subset_size])

  model = ModelStats()
  model.construct("Pipeline Test",
                  data,
                  preprocessor=_cleaner,
                  debug=True)

  model.summary()



if __name__=="__main__":
  main()
