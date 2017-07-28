import re
from sys import argv
import csv
import glob


def extract_results(fn, statistic):
    result_line = r"^sigma=([0-9.]+) omega=([0-9.]+) Average {}: ([0-9.]+)".format(statistic)
    results = []
    with open(fn) as f:
        for line in f:
            res = re.findall(result_line, line)
            if res:
                results.append(map(float, res[0]))

    return results


def extract_and_write(fn, csv_fn, statistic):
    results = extract_results(fn, statistic)

    sigmas = sorted(set([res[0] for res in results]))
    omegas = sorted(set([res[0] for res in results]))

    results_dict = {sigma: {} for sigma in sigmas}

    for res in results:
        sigma, omega, stat = res
        results_dict[sigma][omega] = stat

    with open(csv_fn, 'w') as f:
        writer = csv.DictWriter(f, ['sigma'] + omegas)
        writer.writeheader()
        for sigma in sigmas:
            row = {'sigma': sigma}
            row.update(results_dict[sigma])
            writer.writerow(row)

if __name__ == '__main__':
    try:
        statistic = argv[1]
    except IndexError:
        statistic = 'BCR'

    files = glob.glob('PUK_SCAN*.out')

    for fn in files:
        csv_fn = fn.replace('.out', '.csv')

        extract_and_write(fn, csv_fn, statistic)
