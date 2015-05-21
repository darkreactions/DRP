import csv, json,sys,os


def process_cg(cg_path, output_path="restart.json", sdf_path="sdf", jchem_path = "/home/praccugl/ChemAxon/JChem/bin"):
    import calculate_cg_entry
    CGEntry = calculate_cg_entry.CGEntry

    count = 0
    cg_properties = {}

    if not os.path.exists(sdf_path):
        os.makedirs(sdf_path)

    with open(cg_path) as cg:
	cg.readline() # waste the title
        for line in csv.reader(cg):
            sim, smiles, abbrev = line
            m_type = "Org"
            if abbrev in cg_properties:
                print ("Repeated abbrev! OHNOES: {0}".format(abbrev))
                continue
            if smiles == "":
                print "No smiles for {0}".format(abbrev)
                continue
            fake_abbrev = str(count) + filter(str.isalnum, abbrev)
            try:
                cg_properties[abbrev] = CGEntry(abbrev, fake_abbrev, smiles, m_type, jchem_path,sdf_path).get_properties() 
                cg_properties[abbrev]["sim"] = sim
            except Exception as e:
                print str(e) + " failure in " + abbrev
            count += 1
    json.dump(cg_properties, open(output_path,"w"))
            

def main():
    args = ["cg.csv", "restart.json", "sdf", "/home/praccugl/ChemAxon/JChem/bin"]
    if len(sys.argv) > 1:
        for idx in range(len(sys.argv) - 1):
            args[idx] = sys.argv[idx+1]
    print args
    process_cg(args[0], args[1], args[2], args[3])
            

if __name__ == "__main__":
    main()
