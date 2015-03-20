import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

def get_compound_map():
  def clean_row(row):
    cleaned = []
    for elem in row:
      if elem and "\""==elem[0] and "\""==elem[-1]: elem = elem[1:-1]
      cleaned.append(elem)
    return cleaned

  import csv

  with open("raw/amines.csv") as f:
    matrix = [clean_row(line) for line in csv.reader(f)]

  headers = matrix.pop(0)
  comp_header = "CHEMICAL_NAME"
  comp_index = headers.index(comp_header)

  return {row[comp_index]:{header:row[i] for i, header in enumerate(headers)}
                                       for row in matrix}


def update_compounds():
  from DRP.models import Recommendation
  compound_map = get_compound_map()

  recs = Recommendation.objects.filter()



if __name__=="__main__":
  update_compounds()
