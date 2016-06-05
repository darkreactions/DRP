import os
import sys
full_path = os.path.dirname(os.path.realpath(__file__)) + "/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
    sys.path = [django_path] + sys.path
    os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


DEBUG = True
PRINT_DETAILS = False

# Seeds: ms115.6 jho213.20 jho148.2 jho252.5 (NOTE)
#SPECIFIC_SEED = "ms115.6"
#SEED = SPECIFIC_SEED.split(".")[0]
SEED = "Spreadsheet"
SUFFIX = ""

# If `universe is 'None', then the default universe is specified in the metric.
universe = None

# A metric that returns 1 if two points are in the same "category" or else 0.


def same_category_metric(point, other):
    if "Se" in point.atoms and "Se" in other.atoms:
        return 1
    elif "Te" in point.atoms and "Te" in other.atoms:
        return 1
    else:
        return 0

# Prepare the metric.
from DRP.recommendation.metrics import get_default_metric
metric = get_default_metric(universe=universe, debug=DEBUG)


from DRP.model_building.rxn_calculator import headers

# The field to use to make sure a point isn't compared to itself.
id_field = "XXXtitle"  # TODO: This should change to "ref", but low-priority.
id_index = headers.index(id_field)


distance_cache = {}
calc_cache = {}


def distance(point, other):

    if point not in calc_cache:
        calc_cache[point] = point.get_calculations_list()
    point_calcs = calc_cache[point]

    if other not in calc_cache:
        calc_cache[other] = other.get_calculations_list()
    other_calcs = calc_cache[other]

    key = (point_calcs[id_index], other_calcs[id_index])
    reverse_key = (other_calcs[id_index], point_calcs[id_index])

    if reverse_key in distance_cache:
        return distance_cache[reverse_key]

    else:

        if not key in distance_cache:
            distance_cache[key] = metric(point_calcs, other_calcs)

        return distance_cache[key]


def filter_out_identical(point, others):
    from DRP.models import Data
    if isinstance(point, Data):
        others = filter(lambda other: other.ref != point.ref, others)
    else:
        point_id = point[id_index]
        others = filter(lambda other: other[id_index] != point_id, others)
    return others


def get_knn_tuples(point, others, k):
    """
    Gather the k closest (point, distance) tuples.
    """

    others = filter_out_identical(point, others)

    neighbors = [(other, distance(point, other)) for other in others]
    neighbors.sort(key=lambda neighbor: neighbor[1])
    return neighbors[:k]


def exact_knn_distance(point, others, k):
    """
    Return the exact distance from a `point` to its Kth nearest neighbor.
    """

    knn = get_knn_tuples(point, others, k)
    _, knn_dist = knn[k - 1]
    return knn_dist


def average_knn_distance(point, others, k):
    """
    Calculate the average distance to the k-nearest neighbors for some `point`.
    """

    knn = get_knn_tuples(point, others, k)
    distances = map(lambda tup: tup[1], knn)
    total = float(sum(distances))

    return total / len(knn)


def category_knn_distance(point, others, k):

    knn = get_knn_tuples(point, others, k)
    matches = sum([same_category_metric(point, other) for other, dist in knn])
    return matches


def get_research_points():

    # Used to grab the seed spawn.
    from DRP.research.casey.retrievalFunctions import get_data_from_ref_file
    data = get_data_from_ref_file("DRP/research/casey/raw/033115_datums.txt")
    #data = filter(lambda d: "Te" in d.atoms and "V" in d.atoms, data)

    """
  # Used to grab the data for a specific seed's spawn.
  from DRP.research.casey.retrievalFunctions import get_data_from_ref_file
  data = get_data_from_ref_file("DRP/research/casey/raw/{}_spawn.txt".format(SEED))
  """

    """
  # Used to grab the data used in Alex's 03-09-15 spreadsheet.
  from DRP.research.casey.retrievalFunctions import get_data_from_ref_file
  data = get_data_from_ref_file("DRP/research/casey/raw/030915_datums.txt")
  """

    """
  # Used for the seed KNN graphs.
  from DRP.research.casey.retrievalFunctions import get_data_from_ref_file
  data = get_data_from_ref_file("DRP/research/casey/raw/030915_seeds.txt")
  data = filter(lambda d: "Se" in d.atoms and "V" in d.atoms, data)
  """

    """
  # Used for the seed model.
  from DRP.research.casey.retrievalFunctions import get_data_from_ref_file
  data = get_data_from_ref_file("DRP/research/casey/raw/033115_model.txt")
  data = filter(lambda d: "Se" in d.atoms and "V" in d.atoms, data)
  """

    """
  # Used for the seed model.
  from DRP.research.casey.retrievalFunctions import get_data_from_ref_file
  data = get_data_from_ref_file("DRP/research/casey/raw/030915_intuition.txt")
  """

    """
  # Used for the average KNN distance calculations.
  from DRP.retrievalFunctions import get_valid_data
  from DRP.retrievalFunctions import filter_by_date
  data = get_valid_data()
  data = filter_by_date(data, "04-02-2014", "before")
  outcomes = {"1","2"}
  data = filter(lambda entry: entry.outcome in outcomes, data)
  data = filter(lambda d: "Se" in d.atoms and "V" in d.atoms, data)
  """

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
    data = filter_by_date(data, "04-02-2014", "before")
    #outcomes = {"1","2"}
    #data = filter(lambda entry: entry.outcome in outcomes, data)

    """
  data = filter(lambda d: "Se" in d.atoms and "V" in d.atoms, data)
  data = [d.get_calculations_list(debug=True) for d in data]
  """

    """
  # Graph a single reaction with a specific ref.
  data = get_valid_data()
  data = filter(lambda datum: datum.ref.lower() == "jho252.5", data)
  data = [d.get_calculations_list(debug=True) for d in data]
  """

    return data


def get_knn_research_results(k_range, mode):

    if DEBUG:
        print "Gathering research points..."
    points = get_research_points()
    if not points:
        raise Exception("No research points found!")

    if DEBUG:
        print "Gathering other research points..."
    others = get_research_others()
    if not others:
        raise Exception("No \"other\" research points found!")

    results = {k: [] for k in k_range}

    print "Expanding points..."
    for other in others:
        calc_cache[other] = other.get_calculations_list()

    for point in points:
        calc_cache[point] = point.get_calculations_list()

    if mode == "exact":
        knn_distance = exact_knn_distance

    elif mode == "average":
        knn_distance = average_knn_distance

    elif mode == "category-composition":
        knn_distance = category_knn_distance

    else:
        raise Exception("Unknown `mode` specified!")

    for k in k_range:
        print "{} k={} distance...".format(mode, k)
        for i, point in enumerate(points):
            result = knn_distance(point, others, k)
            results[k].append((point, result))

    return results


def calculate_avg_distance(low, high):
    from DRP.graph import get_graph

    k_range = xrange(low, high + 1)
    results = get_knn_research_results(k_range, "exact")

    for k, reactions in results.items():
        dists = sorted([dist for p, dist in reactions])

        # Graph Options
        padding = 0.01  # percent of graph to use as padding.
        num_major_ticks = 10.0
        num_minor_ticks = 50.0

        max_dist = max(dists)
        min_dist = min(dists)

        # Calculate padding for the graph.
        pre_tick_dist = (max_dist - min_dist) / num_major_ticks
        top = max_dist * (1 + pre_tick_dist * padding)
        bottom = min_dist * (1 - pre_tick_dist * padding)
        if bottom < 0:
            bottom = 0

        line = {"K={}".format(k): dists}
        percentage_baseline = range(len(dists))

        graph = get_graph(line, percentage_baseline,
                          xLabel="# of Reactions",
                          yLabel="Exact KNN Distance",
                          tick_range=(bottom, top),
                          major_tick=(top - bottom) / num_major_ticks,
                          minor_tick=(top - bottom) / num_minor_ticks,
                          show_legend=True,
                          show_minor=True,
                          show_mean=True,
                          )

        graph.show()
        raw_input("Press Enter to continue.")


def get_k_avgs(csv_filename_or_content, mode):
    import csv

    exact_K = 50
    avg_K = 50

    if isinstance(csv_filename_or_content, str):
        with open(csv_filename_or_content, "r") as f:
            reader = csv.reader(f)
            matrix = [row for row in reader]
            data = matrix[1:]

    else:
        data = csv_filename_or_content[1:]

    if mode == "exact":
        data = [row[7:7 + exact_K] for row in data]

    elif mode == "average":
        data = [row[7 + avg_K:] for row in data]

    else:
        raise Exception("Invalid mode specified!")

    totals = [0.0 for elem in data[0]]
    for row in data:
        for i, elem in enumerate(row):
            totals[i] += float(elem)

    num_entries = len(data)
    avgs = [total / num_entries for total in totals]

    return avgs


def write_bucket_to_CSV(filename, bucket, x_axis, mode):
    import csv

    with open(filename, "w") as f:
        writer = csv.writer(f)

        headers = ["Series"] + map(lambda k: "k={} {}".format(k, mode), x_axis)
        writer.writerow(headers)

        for label, values in bucket.items():
            row = [label] + values
            writer.writerow(row)

    print "KNN chart CSV written: {}".format(filename)


def knn_research_graphs(low, high):
    from DRP.graph import get_graph

    mode = "exact"

    k_range = xrange(low, high + 1)

    prefix = "results/knn_calculations_"

    # Gather the Se and Te reactions.
    import csv

    outcome = 2
    is_te = 3
    is_se = 4

    with open(prefix + "all_1To50.csv") as f:
        reader = csv.reader(f)
        data = [row for row in reader]
        headers = data.pop(0)

        bad_outcomes = {"1", "2"}
        good_outcomes = {"3", "4"}

        te = filter(lambda row: row[is_te] == "True", data)
        se = filter(lambda row: row[is_se] == "True", data)

        te_bad = filter(lambda row: row[outcome] in bad_outcomes, te)
        te_good = filter(lambda row: row[outcome] in good_outcomes, te)

        se_bad = filter(lambda row: row[outcome] in bad_outcomes, se)
        se_good = filter(lambda row: row[outcome] in good_outcomes, se)

    buckets = {
        #"Average": get_k_avgs(prefix + "AllSuccess_1To50.csv", mode),
        "Se Fail to Se Failures": get_k_avgs(prefix + "SeFailuresToSeFailures_1To50.csv", mode),
        "Se Fail to Te Success": get_k_avgs(prefix + "SeFailuresToTeSuccesses_1To50.csv", mode),
        #"Te Fail": get_k_avgs([headers]+te_bad, mode),
        #"Te Success": get_k_avgs([headers]+te_good, mode),
        #"Se Fail": get_k_avgs([headers]+se_bad, mode),
        #"Se Success": get_k_avgs([headers]+se_good, mode),
    }

    """
  with open(prefix + "SeSuccess_1To50.csv") as f:
    reader = csv.reader(f)
    data = [row for row in reader]
    headers = data.pop(0)
    ref = 1
    for row in data:
      buckets[row[ref]] = get_k_avgs([headers]+[row], mode)

  with open(prefix + "TeSuccess_1To50.csv") as f:
    reader = csv.reader(f)
    data = [row for row in reader]
    headers = data.pop(0)
    ref = 1
    for row in data:
      buckets[row[ref]] = get_k_avgs([headers]+[row], mode)
  """

    # Graph Options
    padding = 0.01  # percent of graph to use as padding.
    num_major_ticks = 10.0
    num_minor_ticks = 50.0

    max_dist = 0.0
    min_dist = float("inf")
    for point, dists in buckets.items():
        for dist in dists:
            if dist > max_dist:
                max_dist = dist
            if dist < min_dist:
                min_dist = dist

    # Calculate padding for the graph.
    pre_tick_dist = (max_dist - min_dist) / num_major_ticks
    top = max_dist * (1 + pre_tick_dist * padding)
    bottom = min_dist * (1 - pre_tick_dist * padding)
    if bottom < 0:
        bottom = 0

    x_range = list(k_range)
    write_bucket_to_CSV("results/KNN_distance_chart_{}.csv".format(mode), buckets, x_range, mode)

    graph = get_graph(buckets, x_range,
                      xLabel="# Nearest Neighbors (K)",
                      yLabel="{} Distance of K Nearest Neighbors".format(mode.capitalize()),
                      tick_range=(bottom, top),
                      major_tick=(top - bottom) / num_major_ticks,
                      minor_tick=(top - bottom) / num_minor_ticks,
                      show_legend=True,
                      show_minor=True
                      )
    graph.show()
    raw_input("Press Enter to continue...")


def matrix_to_csv(matrix, filename):
    import csv
    with open(filename, "w") as f:
        writer = csv.writer(f)
        for line in matrix:
            writer.writerow(line)

    print "'{}' written!".format(filename)


def make_distance_csv(low, high):
    k_range = xrange(low, high + 1)
    exact_results = get_knn_research_results(k_range, "exact")
    avg_results = get_knn_research_results(k_range, "category-composition")

    # Load the sets of points so we can check which point belongs in which class.
    with open(django_path + "/DRP/research/casey/raw/030915_intuition.txt") as f:
        int_set = f.read().lower().replace("\n", "").split(" ")
    with open(django_path + "/DRP/research/casey/raw/033115_model.txt") as f:
        model_set = f.read().lower().replace("\n", "").split(" ")

    # Model seed sets.
    with open(django_path + "/DRP/research/casey/raw/ms115_spawn.txt") as f:
        ms115_set = f.read().lower().replace("\n", "").split(" ")
    with open(django_path + "/DRP/research/casey/raw/jho148_spawn.txt") as f:
        jho148_set = f.read().lower().replace("\n", "").split(" ")
    with open(django_path + "/DRP/research/casey/raw/jho213_spawn.txt") as f:
        jho213_set = f.read().lower().replace("\n", "").split(" ")
    with open(django_path + "/DRP/research/casey/raw/jho252_spawn.txt") as f:
        jho252_set = f.read().lower().replace("\n", "").split(" ")

    final = {}
    for k, point_tups in exact_results.items():
        avg_tups = avg_results[k]
        for (point, exact_dist), (point2, avg_dist) in zip(point_tups, avg_tups):
            if point != point2:
                raise Exception("Order is wrong!")

            # Get the ref of the seed from which this reaction was spawned.
            ref = point.ref.lower()
            seed = ""
            if ref in ms115_set:
                seed = "ms115.6"
            if ref in jho148_set:
                seed = "jho148.2"
            if ref in jho213_set:
                seed = "jho213.20"
            if ref in jho252_set:
                seed = "jho252.5"

            if point not in final:
                final[point] = {
                    "outcome": point.outcome,
                    "ref": point.ref,
                    "seed": seed,
                    "Te": "Te" in point.atoms,
                    "Se": "Se" in point.atoms,
                    "Intuition": ref in int_set,
                    "Model": ref in model_set,
                }

            final[point]["Exact Distance K={}".format(k)] = exact_dist
            final[point]["Category Composition K={}".format(k)] = avg_dist

    keys = [key for key in final[final.keys()[0]].keys() if "Dist" not in key]
    columns = sorted(keys, reverse=True)
    columns += ["Exact Distance K={}".format(k) for k in k_range]
    columns += ["Category Composition K={}".format(k) for k in k_range]

    matrix = [columns]
    matrix += [[calcs[col] for col in columns] for point, calcs in final.items()]

    filename = "knn_calculations_{}_{}To{}{}.csv".format(SEED, low, high, SUFFIX)
    filepath = "{}/DRP/research/casey/results/{}".format(django_path, filename)
    matrix_to_csv(matrix, filepath)


def get_research_point_time_range():
    data = get_research_points()
    data.sort(key=lambda datum: datum.creation_time_dt)
    print list(data)[0].creation_time_dt
    print list(data)[-1].creation_time_dt


def main():
    #knn_research_graphs(1, 50)
    # calculate_avg_distance(1,25)
    make_distance_csv(1, 1)

if __name__ == "__main__":
    main()
