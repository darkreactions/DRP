import sys, csv, os, random

def help_info():
        print "You need arguments! Usage:\n  Regular use: script.py <source_file> <new_prefix>"
        print "If you want to disable shuffle, use: script.py <source_file> <new_prefix> 1\n  For this help: script.py help"

def remove_col(row, outcome):
    if outcome:
        return row[:-2] + ["1" if row[-1] not in ['1','2','3','4',1,2,3,4] else row[-1]]
    else:
        return row[:-2] + ["1" if row[-2] not in ['1','2',1,2] else row[-2]]

def preface(headers, outcome, prefix, specials):
    res = "%  COMMENT \n%  NAME, DATE\n@relation " + prefix 
    for header in headers:
        if header in specials.keys():
            if not (outcome and header == "purity") and not (not outcome and header == "outcome"):
                res += "\n@ATTRIBUTE " + header + " " + specials[header]
        else:
            res += "\n@ATTRIBUTE " + header + " NUMERIC"
    res += "\n\n@DATA\n"
    return res


def process_rows(rows, headers, prefix, outcome=True):
    specials = gen_specials() 
    XXX = 0
    for header in headers:
        if header[0:3] == "XXX":
            XXX += 1
        headers = headers[XXX:]
    result = preface(headers, outcome, prefix, specials) 
    for row in rows:
        result +=  ",".join(remove_col(row[XXX:], True)) + "\n"
    return result
        

def gen_specials():
    specials = {"outcome": "{1,2,3,4}", "slowCool": "{yes,no}", 
            "leak": "{yes,no}", "purity": "{1,2}"}
    import rxn_calculator
    for atom in rxn_calculator.atomsz + rxn_calculator.bools:
        specials[atom] = "{yes,no}"
    return specials

def clean(prefix, specials = gen_specials()):
    XXX = 0
    headers = None 
    do_shuffle = True
    if len(sys.argv) == 3:
        do_shuffle = False
    with open(prefix+".csv","rb") as r:
        reader = csv.reader(r)
        headers = reader.next()
        for header in headers:
            if header[0:3] == "XXX":
                 XXX += 1
        headers = headers[XXX:]
        rows = [ row for row in reader]
        random.shuffle(rows)
        reader = rows
        with open(prefix + "_out.arff","w") as out:
            out.write(preface(headers, True, prefix, specials))
            for row in reader:
                out.write(",".join(remove_col(row[XXX:], True)) + "\n")

if __name__ == "__main__":
    if sys.argv[1] in ["help","h","-h","--help","--h","-help"]:
        help_info()
    elif False: # TODO: Filter out bad sys.argv[1], sys.argv[2], and make sure the files aren't already touched, and allow notouch to be overridden
        pass
    else:
        specials = gen_specials()
        
        #CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        #PREFIX_DIR = "{0}/{1}/".format(CHEMML_DIR, sys.argv[1])    
        #prefix = "{0}/{1}/{1}".format(CHEMML_DIR, sys.argv[1])    
        prefix = sys.argv[1]
        clean(prefix, specials)
