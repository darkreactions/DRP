import glob
from sys import argv
from numpy import mean

def get_val_from_file(filepath, name):
    with open(filepath, 'rb') as f:
        raw_lines = f.readlines()
    return get_val(raw_lines, name)
            
def get_val(lines, name):
    start_line = name + ':'
    for line in lines:
        if line.startswith(start_line):
            val = line.split()[1]
            return val

def get_average_val_from_file(filepath, name):
    with open(filepath, 'rb') as f:
        raw_lines = f.readlines()
    
    return get_average_val(raw_lines, name)

def get_average_val(lines, name):
    start_line = name + ':'
    skipped_first = False
    vals = []
    for line in lines:
        if line.startswith(start_line):
            if skipped_first:
                val = float(line.split()[1])
                vals.append(val)
            else:
                skipped_first = True
    return mean(vals)

if __name__ == '__main__':  
    directory = 'legacy_tests/model_out/'

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
        else:
            raise RuntimeError("Don't know how to handle this many parts in a filenmae")
        
        
        print "{}\t{}\t{}".format(mod_string, get_val_from_file(fp, 'BCR'), get_val_from_file(fp, 'Accuracy'))
