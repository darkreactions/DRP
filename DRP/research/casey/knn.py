import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


debug = True


#If `universe is 'None', then the default universe is specified in the metric.
universe = None

# Prepare the metric.
from DRP.recommendation.metrics import get_default_metric
distance = get_default_metric(universe=universe, debug=debug)


def filter_out_identical(point, others):

  from DRP.model_building.rxn_calculator import headers

  # The field to use to make sure a point isn't compared to itself.
  id_field = "ref"
  id_index = headers.index(id_field)

  point_id = point[id_index]
  others = filter(lambda other: other[id_index]!=point_id, others)

  return others


def get_knn_tuples(point, others, k):
  """
  Gather the k closest (point, distance) tuples.
  """

  others = filter_out_identical(point, others)

  neighbors = [(other, distance(point,other)) for other in others]
  neighbors.sort(key=lambda neighbor: neighbor[1])
  return neighbors[:k]


def average_knn_distance(point, others, k):
  """
  Calculate the average distance to the k-nearest neighbors for some `point`.
  """

  knn = get_knn_tuples(point, others, k)
  distances = map(lambda tup: tup[1], knn)
  total = float(sum(distances))
  return total/len(knn)


def get_research_points():
  from DRP.retrievalFunctions import get_valid_data

  data = list(get_valid_data())[:5]

  return data


def get_research_others():
  from DRP.retrievalFunctions import get_valid_data

  data = [p.get_calculations_list() for p in get_valid_data()[:100]]

  return data


def main():

  k = 1

  if debug: print "Gathering research points..."

  points = get_research_points()
  others = get_research_others()

  print "Average k={} distance...".format(k)
  for i, point in enumerate(points):
    avg_dist = average_knn_distance(point.get_calculations_list(), others, k)
    print "({}) {} : {}".format(i, point.ref, avg_dist)

  if debug: print "-- Finished!"


if __name__=="__main__":
  main()

