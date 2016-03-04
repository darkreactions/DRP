import glob
from sys import argv

def get_val(filepath, name):
    start_line = name + ':'
    with open(filepath, 'rb') as f:
        raw_lines = f.readlines()
    for line in raw_lines:
        if line.startswith(start_line):
            val = line.split()[1]
            return val


if __name__ == '__main__':
    directory = 'model_for_full_draft_figs'
    name = 'Accuracy'

    filepaths = sorted(glob.glob('{}/*'.format(directory)))

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

        
        
        print "{}\t{}".format(mod_string, get_val(fp, name))
