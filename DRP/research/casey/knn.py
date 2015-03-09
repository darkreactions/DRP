import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def distance(point, other):
  from DRP.recommendation.metrics import default_metric
  return default_metric(point, other)


def get_knn_tuples(point, universe, k):
  """
  Gather the k closest (point, distance) tuples.
  """

  neighbors = [(other, distance(point,other)) for other in universe ]
  neighbors.sort(key=lambda neighbor: neighbor[1])
  return neighbors[:k]


def average_knn_distance(point, universe, k):
  """
  Calculate the average distance to the k-nearest neighbors for some `point`.
  """

  knn = get_knn_tuples(point, universe, k)
  distances = map(lambda tup: tup[1], knn)
  total = float(sum(distances))
  return total/len(knn)


def main():
  print "Here!"



if __name__=="__main__":
  main()

