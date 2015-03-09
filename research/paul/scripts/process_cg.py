import csv, json,sys,os


def single_line(line):
    import calculate_cg_entry
    CGEntry = calculate_cg_entry.CGEntry
    print line
    abbrev, compound, cas, m_type, smiles, count, jchem_path, sdf_path = line
    if smiles == "":
        print "No smiles for {0}".format(abbrev)
        return None
    fake_abbrev = str(count) + filter(str.isalnum, abbrev)
    return (abbrev, CGEntry(abbrev, fake_abbrev, smiles, m_type, jchem_path,sdf_path).get_properties())

	

def process_cg(cg_path, output_path="restart.json", sdf_path="sdf", jchem_path = "/home/praccugl/ChemAxon/JChem/bin"):
    if not os.path.exists(sdf_path):
        os.makedirs(sdf_path)


    lines = []
    with open(cg_path) as cg:
        cg.readline() # waste the title
        reader = csv.reader(cg)
        json.dump(process_cg_rows(reader), open("restart.json","w"))

def process_cg_rows(reader):
    count = 0
    cg_properties = {}
    cnt  = 0
    for l in reader:
        lines.append(l + [str(cnt), jchem_path, sdf_path ])
    import multiprocessing
    pool = multiprocessing.Pool(processes=10)
    results = pool.map(single_line, lines)
    cg_properties = {r[0]:r[1] for r in results if r is not None}
    return cg_properties

            

def main():
    args = ["cg.csv", "restart.json", "sdf", "/home/praccugl/ChemAxon/JChem/bin"]
    if len(sys.argv) > 1:
        for idx in range(len(sys.argv) - 1):
            args[idx] = sys.argv[idx+1]
    process_cg(args[0], args[1], args[2], args[3])
            

if __name__ == "__main__":
    import time
    print time.time()
    main()
    print time.time()
