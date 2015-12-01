from DRP.settings import BASE_DIR

import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

def default_preprocessor(data):
  # Takes a list of Data objects prepended by the headers
  headers = data.pop(0)

  expanded = []
  for d in data:
    try:
      expanded.append(d.get_calculations_list())
    except:
      pass

  data = [headers] + expanded

  return default_calc_list_preprocessor(data)

def default_row_preprocessor(data):
  # EXPERIMENTAL, not yet working
  # Takes a list of raw rows of the features prior to expansion, prepended by the headers

  from DRP.model_building import parse_rxn, load_cg
  import json
  cg_props = load_cg.get_cg()
  ml_convert = json.load(open("{}/DRP/model_building/mlConvert.json".format(BASE_DIR)))
  #ml_convert = json.load(open(django_path+"/DRP/model_building/mlConvert.json"))

  headers = data.pop(0)

  expanded = []
  for d in data:
    try:
      d_expanded = parse_rxn.parse_rxn(d, cg_props, ml_convert)
      expanded.append(d_expanded)
    except:
      pass

  data = [headers] + expanded

  return default_calc_list_preprocessor(data)

def default_calc_list_preprocessor(data):
  # Takes a list of calculations lists prepended by the headers
  # (this is the internals of the default_preprocessor, where all the work is done)
  def is_num(elem):
    try:
      float(elem)
      return True
    except:
      return False

  def norm(elems):
    import operator
    elems = map(float,elems)
    avg = float(reduce(operator.add, elems))/len(elems)

    return [elem/avg for elem in elems]

  # Variable Setup
  normalize = False

  headers = data.pop(0)

  # Set outcomes of 0 to be 1.
  outcome_index = headers.index("outcome")
  outcomes = [int(row[outcome_index]) for row in data]
  outcomes = [elem if elem>0 else 1 for elem in outcomes]
  for i, outcome in enumerate(outcomes):
    data[i][outcome_index] = str(outcome)


  # Remove weird unicode characters.
  for i, entry in enumerate(data):
    for j, elem in enumerate(entry):
      if type(elem)==str or type(elem)==unicode:
        data[i][j] = elem.replace(u"\u2019", "'")


  if normalize:

    # Normalize everything except the categorized columns.
    categories = [
        "outcome",
    ]
    cat_indexes = {headers.index(cat) for cat in categories if cat in headers}

    cols = [[unicode(row[i]) for row in data] for i, h in enumerate(headers)]
    cols = [norm(col) if (all(map(is_num, col)) and i not in cat_indexes) else col
                    for i, col in enumerate(cols)]
    data = [[col[i] for col in cols] for i in xrange(len(data))]

  return [headers]+data


# TODO: Un-tested and perhaps unnecessary...
def purging_preprocessor(data):
  def remover(row):
    to_remove = {
      "water"
    }
    return not any(unicode(elem).lower() in to_remove for elem in row)

  headers = data.pop(0)
  data = [d.get_calculations_list() for d in data]

  data = filter(remover, data)

  return [headers] + data


def category_preprocessor(data):
  def categorize(data):
    return "Te" in data.atoms and "V" in data.atoms

  new_data = default_preprocessor(data)

  headers = new_data.pop(0)
  new_data = [row+[categorize(datum)] for row, datum in zip(new_data, data)]

  return [headers] + new_data

