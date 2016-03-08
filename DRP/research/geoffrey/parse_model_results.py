import glob
from sys import argv
from numpy import mean

def get_val(filepath, name):
    start_line = name + ':'
    with open(filepath, 'rb') as f:
        raw_lines = f.readlines()
    for line in raw_lines:
        if line.startswith(start_line):
            val = line.split()[1]
            return val

def get_average_val(filepath, name):
    start_line = name + ':'
    with open(filepath, 'rb') as f:
        raw_lines = f.readlines()
    skipped_first = False
    vals = []
    for line in raw_lines:
        if line.startswith(start_line):
            if skipped_first:
                val = float(line.split()[1])
                vals.append(val)
            else:
                skipped_first = True
    return mean(vals)


if __name__ == '__main__':
    directory = 'legacy_tests'

    filepaths = sorted(glob.glob('{}/*.out'.format(directory)))


    print "\t\t\t{}\t{}".format('BCR', 'Accuracy')
    for fp in filepaths:
        fn = fp.split('/')[-1]
        fn_split = fn[:-4].split('_')

        if len(fn_split) == 1:
            model = fn_split[0]
            mod_string = 'none\tEuclidean\t{}'.format(model)
        elif len(fn_split) == 2:
            model = fn_split[0]
            fs = fn_split[1]
            mod_string = '{}\tEuclidean\t{}'.format(fs, model)
        elif len(fn_split) == 3:
            model = fn_split[0]
            metric = fn_split[1]
            fs = fn_split[2]
            mod_string = '{}\t{}\t{}'.format(fs, metric, model)

        
        
        print "{}\t{}\t{}".format(mod_string, get_average_val(fp, 'BCR'), get_average_val(fp, 'Accuracy'))
