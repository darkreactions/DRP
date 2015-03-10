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
metric = get_default_metric(universe=universe, debug=debug)


from DRP.model_building.rxn_calculator import headers

# The field to use to make sure a point isn't compared to itself.
id_field = "XXXtitle" #TODO: This should change to "ref", but low-priority.
id_index = headers.index(id_field)


distance_cache = {}
def distance(point, other):

  # Cache the distance between points for future use.
  key = (point[id_index], other[id_index])
  reverse_key = (other[id_index], point[id_index])

  if reverse_key in distance_cache:
    return distance_cache[reverse_key]

  else:

    if not key in distance_cache:
      distance_cache[key] = metric(point, other)

    return distance_cache[key]



def filter_out_identical(point, others):
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

  if point[0] in {"jk252.5","jho252.5"} :
    print "{} Nearest Neighbors of {}:".format(k, point[0])
    print "\n".join(["\t{} : {}".format(tup[0][0], tup[1]) for tup in knn])

  return total/len(knn)


def get_research_points():
  from DRP.research.casey.retrievalFunctions import get_data_from_ref_file

  data = get_data_from_ref_file("DRP/research/casey/raw/030915_seeds.txt")

  """
  #TODO: quickie-info
  dates = [d.creation_time_dt for d in data]
  dates.sort()
  print "Latest: {}".format(dates[-1])
  print "Earliest: {}".format(dates[0])

  print "Total: {}".format(len(data))
  te = filter(lambda d: "Te" in d.atoms, data)
  se = filter(lambda d: "Se" in d.atoms, data)
  neither = filter(lambda d: not "Te" in d.atoms and not "Se" in d.atoms, data)
  print "Te: {}".format(len(te))
  print "Se: {}".format(len(se))
  print "neither: {}".format(len(neither))
  """

  return data


def get_research_others():
  from DRP.retrievalFunctions import get_valid_data
  from DRP.retrievalFunctions import filter_by_date

  data = get_valid_data()
  data = filter_by_date(data, "05-21-2014", "before")
  data = [d.get_calculations_list(debug=True) for d in data]

  return data



def knn_research_graphs(low, high):
  from DRP.graph import get_graph

  if debug: print "Gathering research points..."
  points = get_research_points()
  if not points: raise Exception("No research points found!")

  if debug: print "Gathering other research points..."
  others = get_research_others()
  if not others: raise Exception("No \"other\" research points found!")

  k_range = xrange(low, high+1)

  results = {k:[] for k in k_range}
  calc_cache = {}


  for k in k_range:
    print "Average k={} distance...".format(k)
    for i, point in enumerate(points):

        # Store the `calculations_list` of each point for speed-up.
        if point not in calc_cache:
          calc_cache[point] = point.get_calculations_list()

        avg_dist = average_knn_distance(calc_cache[point], others, k)
        results[k].append( (point, avg_dist) )

  # Sort the reactions and their distances into Se/Te buckets.
  buckets = {"Te":{}, "Se":{}}

  for k, reactions in results.items():
    for point, dist in reactions:
      if "Te" in point.atoms:
        key = "Te"
      elif "Se" in point.atoms:
        key = "Se"

      if point not in buckets[key]:
        buckets[key][point] = [None for i in k_range]

      buckets[key][point][k-1] = dist


  for key, bucket in buckets.items():
    if bucket:
      print "Graphing {}... ({})".format(key, len(bucket))

      max_dist = 0
      for point, dists in bucket.items():
        for dist in dists:
          if dist>max_dist: max_dist = dist

      graph = get_graph(bucket, list(k_range),
                        xLabel="# Nearest Neighbors (K)",
                        yLabel="Average KNN Distance",
                        tick_range=(0,max_dist),
                        major_tick=max_dist/10.0,
                        show_legend=True,
                        show_minor=False)
      graph.show()
      raw_input("Press enter to continue...")
    else:
      print "Skipping empty bucket `{}`...".format(key)


def main():
  knn_research_graphs(1,10)



if __name__=="__main__":
  main()

