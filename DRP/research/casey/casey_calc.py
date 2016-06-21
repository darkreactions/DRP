#!/usr/bin/python

import sys

# give 'all' to count all se and te or all outcomes


def countif(all_data, lower_dist, upper_dist, k, avg_exact, se_te, desired_outcome, model_intuition):
    count = 0
    for data in all_data:
        outcome = int(data[2])
        isTe = data[3]
        isSe = data[4]
        fromModel = data[5]
        exactDists = data[7:7 + exact_K]
        avgDists = data[7 + avg_K:]

        if k > exact_K or k > avg_K:
            raise Exception("Data for K={} not found in CSV!".format(k))

        if model_intuition == 'model' and fromModel != 'True':
            continue
        elif model_intuition == 'intuition' and fromModel != 'False':
            continue

        if avg_exact == 'avg':
            dist = float(avgDists[k - 1])
        elif avg_exact == 'exact':
            dist = float(exactDists[k - 1])
        else:
            print "avg_exact must be avg or exact"
            return

        if dist < lower_dist or dist >= upper_dist:
            continue

        if se_te == 'Se' and not isSe == 'True':
            continue
        if se_te == 'Te' and not isTe == 'True':
            continue
        # passing in 'all' or anything else will count everything

        if desired_outcome != 'all' and outcome != desired_outcome:
            continue

        count = count + 1

    return count


def experiment_equal_sized_buckets(all_data, min, max, bucket_size, k, se_te, model_intuition, rows):
    while min < max:
        upper = min + bucket_size
        total_for_range = countif(
            all_data, min, upper, k, 'avg', se_te, 'all', model_intuition)
        rows[0].append(str(min) + '-' + str(upper))
        for outcome in range(1, 5):
            total_for_outcome = countif(
                all_data, min, min + bucket_size, k, 'avg', se_te, outcome, model_intuition)
            if total_for_outcome == 0 or total_for_range == 0:
                percent_for_outcome = 0
            else:
                percent_for_outcome = (
                    float(total_for_outcome) / float(total_for_range)) * 100
            rows[outcome].append(str(percent_for_outcome))
        min = min + bucket_size


def experiment_equal_num_in_buckets(all_data, min, max, num_buckets, k, se_te, model_intuition):
    import math

    row0 = ['k=' + str(k)]
    row1 = ['1']
    row2 = ['2']
    row3 = ['3']
    row4 = ['4']
    rows = []
    rows.append(row0)
    rows.append(row1)
    rows.append(row2)
    rows.append(row3)
    rows.append(row4)

    # Collect Bucket Sizes (by avg KNN distance)
    avg_k_index = 7 + exact_K + (k - 1)
    dists = sorted([float(row[avg_k_index]) for row in all_data])
    bucket_size = int(math.ceil(float(len(dists) / (num_buckets - 1))))
    bucket_sizes = [dists[i] for i in xrange(
        bucket_size, len(dists), bucket_size)] + [dists[-1]]

    if dists[0] > min:
        min = dists[0]

    trimmer = 6
    ignore_boundaries = True

    while (min < max or ignore_boundaries) and bucket_sizes:
        upper = bucket_sizes.pop(0)

        rows[0].append(str(min)[:trimmer] + '-' + str(upper)[:trimmer])
        total_for_range = countif(
            all_data, min, upper, k, 'avg', se_te, 'all', model_intuition)

        for outcome in range(1, 5):
            total_for_outcome = countif(
                all_data, min, upper, k, 'avg', se_te, outcome, model_intuition)
            if total_for_outcome == 0 or total_for_range == 0:
                percent_for_outcome = 0
            else:
                percent_for_outcome = (
                    float(total_for_outcome) / float(total_for_range)) * 100
            rows[outcome].append(str(percent_for_outcome))

        min = upper

    for i in range(0, 5):
        print rows[i]
    print "\n"

    return rows


def experiment_sete_fixed(all_data, min, max, num_buckets, k, se_te, model_intuition):

    bucket_size = (max - min) / float(num_buckets)

    row0 = ['k=' + str(k)]
    row1 = ['1']
    row2 = ['2']
    row3 = ['3']
    row4 = ['4']
    rows = []
    rows.append(row0)
    rows.append(row1)
    rows.append(row2)
    rows.append(row3)
    rows.append(row4)

    experiment_equal_sized_buckets(
        all_data, min, max, bucket_size, k, se_te, model_intuition, rows)

    for i in range(0, 5):
        print rows[i]
    print "\n"

    return rows


def make_3d_plot(matrix_tups):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import numpy as np

    matrix_tups = [matrix_tups[0]]

    for matrix, title in matrix_tups:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        plt.title(filename.split("/")[-1])

        groups = matrix.pop(0)[1:]
        outcomes = []

        for i, row in enumerate(matrix):
            outcomes.append(row.pop(0))

            xs = np.arange(len(row))
            ys = map(float, row)

            ax.bar(xs, ys, zs=i, zdir='y', alpha=0.5)

        plt.xticks(np.arange(len(groups)), groups)
        plt.yticks(np.arange(len(matrix)), outcomes)

        ax.set_xlabel('Distances')
        ax.set_ylabel('Outcome')
        ax.set_zlabel('% of Total')

        ax.view_init(elev=70, azim=315)  # Set the camera position.

        plt.show()


def experiment(all_data, min, max, num_buckets, k):

    # supports: experiment_sete_fixed or experiment_equal_num_in_buckets
    bucketer = experiment_equal_num_in_buckets

    rows = []

    for option in ["model", "intuition"]:
        for division in ["all", "Se", "Te"]:
            key = "{} {}".format(division, option)

            matrix = bucketer(all_data, min, max,
                              num_buckets, k, division, option)
            rows.append((matrix, key))

    make_3d_plot(rows)


if len(sys.argv) < 2:
    print "Usage: " + sys.argv[0] + " <csv file_name>"
    exit(1)

filename = sys.argv[1]
f = open(filename, 'r')

all_data = []
exact_K = 0
avg_K = 0

# Discard the headers, but count how many K are allowed.
for line in f:
    data = line.split(",")
    if data[0] == 'ref' or data[0] == "seed":
        for header in data:
            if "exact distance k=" in header.lower():
                exact_K += 1
            if "average distance k=" in header.lower():
                avg_K += 1
        continue
    all_data.append(data)

f.close()

experiment(all_data, 0, 300, 5, 20)
