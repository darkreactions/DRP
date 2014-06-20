import sys, os, json, csv, math 
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.model_building.rxn_calculator as rxn_calculator 
import traceback

from DRP.models import DataCalc 


def fixmL(t):
    o = ''
    for i in t:
        if ord(i) == ord('/') or ord(i) > ord('9') or ord(i) < ord('.'):
            return o
        o += i
    return o

def parse_rxn(row, rxn_table, ml_convert):
    compound = ["x","x","x","x","x"]
    mass = ["-1","-1","-1","-1","-1"]
    unit = ["","","","",""]
    (title, compound[0], mass[0], unit[0], compound[1], mass[1], unit[1],
        compound[2], mass[2], unit[2], compound[3], mass[3],
        unit[3], compound[4], mass[4], unit[4], Tmax,
        time,  pH, slowCool, leak, outcome, purity) = row[:23]

    for i in range(5):
        if compound[i] == "water": continue
        if compound[i] == "" or rxn_table[compound[i]]["type"] in ['pH']: # TODO: fix
            compound[i] = 'x'
            mass[i] = '-1'
        elif "mL" in unit[i]:
            mass[i] = float(mass[i])*ml_convert[compound[i]]
        elif mass[i] == '' or mass[i] == '?':
            mass[i] = '-1'
        if compound[i] == "":
            compound[i] = "x"
    mass = [float(m) for m in mass]
    Tmax = float(Tmax)
    time = float(time)

    pH = float(pH)
    if not (pH > 0 and pH <= 14): #TODO: fix to allow 0 and 14
        raise Exception("pH {0} is out of bounds: {1}".format(str(pH), title))

    if not purity:
        purity = 1
    purity = int(purity)

    if not outcome:
        outcome = 1

    comb = False
    if not comb:
        outcome = int(outcome)
    else:
        outcome = int(str(outcome) + str(purity))

    output = [title]
    isWater = -1 
    keepList = [1,1,1,1,1]
    organicList = [0,0,0,0,0]
    inorganicList = [0,0,0,0,0]
    waterList = [0,0,0,0,0]
    oxalateList = [0,0,0,0,0]
    (nInorg, nOrg, nOxlike) = (0,0,0)
    for i in range(5):
        if compound[i].lower() == "water":
            isWater = i 
        elif compound[i] == "" or compound[i].lower() in ["x","?"]:
            keepList[i] = 0
    for i in range(5):
        if keepList[i]:
            if rxn_table[compound[i]]["type"] == "Inorg":
                inorganicList[i] = 1 
                nInorg += 1
            elif rxn_table[compound[i]]["type"] == "Org":
                organicList[i] = 1 
                nOrg += 1
            elif rxn_table[compound[i]]["type"] == "Ox":
                oxalateList[i] = 1 
                nOxlike += 1
            elif compound[i].lower() == "water":
                waterList[i] = 1 
                isWater = i
            else:
                raise Exception("Unknown type: {0}".format(compound[i]))
    if isWater == -1:
        raise Exception("Lacks water")
    pH = math.ceil(float(pH))

    bar = None
    compoundAcc = [-1,-1,-1,-1,-1]
    compoundDon = [-1,-1,-1,-1,-1]
    compoundMoles = [-1,-1,-1,-1,-1]

    orgDetails = []
    inorgDetails = []
    oxlikeDetails = []

    for i in range(5):
        if not keepList[i]:
            continue
        a = rxn_table[compound[i]]["mw"] 
        assert(a != 0)
        compoundMoles[i] = mass[i]/a
        if isWater == i and not compoundMoles[i]:
            print row
            raise Exception("Water moles zero: {0}, {1}".format(mass[i], a))
        r = [compound[i], mass[i], compoundMoles[i]]
        if inorganicList[i]:
            inorgDetails += r
        elif organicList[i]:
            orgDetails += r
        elif oxalateList[i]:
            oxlikeDetails += r
        elif i != isWater:
            print row
            print inorganicList
            print organicList
            print oxalateList
            print i
            raise Exception("this should never happen... not inorg, not org, not oxalate? check the source code...")
    for i in [0,1,2]: #TODO: rewrite this to not be so dumb
        if len(inorgDetails) + 1 < 9:
            inorgDetails += [-1,-1,-1]
        if len(orgDetails) + 1 < 6:
            orgDetails += [-1,-1,-1]
        if len(oxlikeDetails) + 1 < 3:
            oxlikeDetails += [-1,-1,-1] 
    output += inorgDetails + orgDetails + oxlikeDetails
    output += [Tmax, time, slowCool.lower(), pH, leak.lower(), nInorg, nOrg, nOxlike, nInorg+nOrg+nOxlike]
    compoundProperties = [[] for k in range(5)]
    for j in range(5):
        if not keepList[j] or isWater == j:
            compoundProperties[j] = [-1 for k in range(19)]
            continue
        a = compound[j]
        compoundProperties[j] += rxn_table[compound[j]]["NopH"]
        compoundProperties[j] += rxn_table[compound[j]]["projectionArea"]
        bar = rxn_table[compound[j]]["polsurf"]
        compoundProperties[j] += [bar[k] for k in range(9*int(pH) - 9, 9*int(pH))]
        bar = rxn_table[compound[j]]["msacc"]
        compoundAcc[j] = bar[int(pH) -1] # TODO: what is this?
        #if abs(compoundAcc[j]) < 0.01:
        #    compoundAcc[j] = 0.0001 # TODO: is this _really_ valid?
        bar = rxn_table[compound[j]]["msdon"]
        compoundDon[j] = bar[int(pH) - 1] # TODO: same
        #if abs(compoundDon[j]) < 0.01:
        #    compoundDon[i] = 0.0001 # TODO: is this valid?
        compoundProperties[j] += [compoundAcc[j], compoundDon[j]]
        if len(compoundProperties[j]) != 19:
            raise Exception("Wat?" + str(len(compoundProperties[j])))

    if nOrg > 0:
        output += rxn_calculator.distList(organicList, compoundProperties)
    else:
        output += [-1 for i in range(76)]
    if nOxlike > 0:
        output += rxn_calculator.distList(oxalateList, compoundProperties)
    else:
        output += [-1 for i in range(76)]

    rxn_props_calculator = rxn_calculator.PropertiesCalculator(compoundMoles, inorganicList, compoundAcc, compoundDon, organicList, isWater)
    output += [rxn_props_calculator.inorg_water_mole_ratio(),
        rxn_props_calculator.org_water_mole_ratio(),
        rxn_props_calculator.org_water_acceptor_on_inorg_ratio(),
        rxn_props_calculator.org_water_donor_on_inorg_ratio(),
        rxn_props_calculator.inorg_org_mole_ratio(),
        rxn_props_calculator.not_water_mole_ratio(),
        ]
    

    smiles = []
    atoms = []
    for j in range(5):
        if keepList[j] == 1 and isWater != j:
            atoms += rxn_table[compound[j]]["atoms"]
            smiles.append((rxn_table[compound[j]]["smiles"], compoundMoles[j]) )
	else:
		smiles.append(None)
    atoms = set(atoms)
    for atom in rxn_calculator.atomsz:
        if atom in atoms:
           output.append("yes")
        else:
           output.append("no")

    atom_counts = {} 
    for idx in range(len(keepList)):
        if (not keepList[idx]) or (isWater == idx): continue
        if compound[idx] in rxn_table and "atom_count" in rxn_table[compound[idx]]: #TODO: check name
            atoms_info = rxn_table[compound[idx]]["atom_count"]
        else:
            at_list = rxn_calculator.atoms_from_smiles(smiles[idx][0])
            atoms_info = {a: at_list.count(a) for a in at_list} 

        for atom in atoms_info:
            if atom not in atom_counts:
                atom_counts[atom] = atoms_info[atom]*compoundMoles[idx]
            else:
                atom_counts[atom] += atoms_info[atom]*compoundMoles[idx] 
            

    output += rxn_calculator.atomic_properties(atoms, smiles, atom_counts)
    output.append(purity)
    output.append(outcome)
    # Needs work: trying to compile all expanded data into the DataCalc model 
    result = map(lambda s: "{0:.4f}".format(s) if type(s) == float else s, output)
    datum = Data.objects.get(id=ID)
    newDataCalc = DataCalc(contents)
    newDataCalc.save()
    datum.calculations = newDataCalc
    datum.save() 
    if datum.calculations:

    else:
    jsonText = json.dumps(result) 
    newDataCalcObj = DataCalc(contents=jsonText)
    newDataCalcObj.save() 
    return result 
    #return map(lambda s: "{0:.4f}".format(s) if type(s) == float else s,  output)  

def parse_rxns(data, name, cg_props = None, ml_convert = None):
    results = []
    
    if ml_convert is None:
        ml_convert = json.load(open("mlConvert.json"))
    if cg_props is None:
        cg_props = json.load(open("restart.json")) 
    for line in data:
        results.append(parse_rxn(line, cg_props, ml_convert))
    return results

def find_ranges(rows):
    ranges = {'temp':{'min':10000000,'max':0}}
    for row in rows:
        compound = ["x","x","x","x","x"]
        mass = ["-1","-1","-1","-1","-1"]
        unit = ["","","","",""]
        (title, compound[0], mass[0], unit[0], compound[1], mass[1], unit[1],
            compound[2], mass[2], unit[2], compound[3], mass[3],
            unit[3], compound[4], mass[4], unit[4], Tmax,
            time,  pH, slowCool, leak, outcome, purity) = row[:23]
        try:
            for i in range(len(compound)):
                if compound[i] == "" or compound[i] == "x" or unit[i] != "g":
                    continue
                if compound[i] not in ranges:
                    ranges[compound[i]] = {'massmin':1000000, 'massmax':0, 'pHmax': 0, 'pHmin': 14}

                m = float(mass[i])
                if ranges[compound[i]]['massmin'] > m and m > 0:
                    ranges[compound[i]]['massmin'] = m 
                if ranges[compound[i]]['massmax'] < m:
                    ranges[compound[i]]['massmax'] = m
                pH = float(pH)
                if ranges[compound[i]]['pHmax'] < pH:
                    ranges[compound[i]]['pHmax'] = pH
                if ranges[compound[i]]['pHmin'] > pH:
                    ranges[compound[i]]['pHmin'] = pH
        except Exception as e:
            print "bork: %s" % str(e)

    return ranges

def find_combinations(rows):
    combos = set()
    for row in rows:
        compounds = tuple(sorted(filter(lambda x: x != "", [row[1], row[4], row[7], row[10], row[13]])))
        combos.add(compounds)
    return combos
    


def build_prior_rows(mass_ranges, row):
    if len(mass_ranges) == 0:
        raise Exception("WTF")
    postfixes = [[mass_ranges[0][0], mass_ranges[0][1] * .5, 'g'], [mass_ranges[0][0], mass_ranges[0][1] * .25, 'g'],
            [mass_ranges[0][0], mass_ranges[0][2] * 1.5, 'g'],  [ mass_ranges[0][0], mass_ranges[0][2] * 2, 'g'] ]
    rows = [row + postfixes[0], row + postfixes[1], row + postfixes[2], row + postfixes[3]]
    if len(mass_ranges) == 1:
        if (len(rows[0]) - 1) / 3 < 5:
            diff = 5 - (len(rows[0]) - 1) / 3
            for i in range(diff):
                rows = [row + ['','',''] for row in rows]

        return rows
    else:
        mr_short = mass_ranges[1:]
        l = []
        for r in rows:
             l += build_prior_rows(mr_short, r)
        return l
    
        

def build_prior(rows):
    prior = []
    ranges = find_ranges(rows)
    combos = find_combinations(rows)
    print ranges
    for combo in combos:
        try:
            pH_max = 0
            pH_min = 14
            mass_ranges = []
            for c in combo:
                mass_ranges.append((c, ranges[c]['massmin'], ranges[c]['massmax']))
                pH_max = max(pH_max, ranges[c]['pHmax']) if ranges[c]['pHmax'] < 14 else 14
                pH_min = min(pH_min, ranges[c]['pHmin'])
            rows = [ r + [170, 170] for r in build_prior_rows(mass_ranges,['--prior'])]
            rows2 = []
            if pH_max < 14 and pH_max > 0:
                for i in range(pH_max + 1, 14 + 1):
                    rows2 += [r + [i] for r in rows]
            else:
                rows2 += [ r + [max(1,pH_max)] for r in rows]
            if pH_min > 1 and pH_min != 14: 
                for i in range(1, pH_min):
                    rows2 += [r + [i] for r in rows]
            else:
                rows2 += [r  + [max(1,pH_min)] for r in rows]
            rows = [r + ['Yes','No',0,1] for r in rows2]
        except Exception as e:
            print str(e)
        prior += rows
    return prior


def main():
    if len(sys.argv) != 4:
        print "not enough args"
        return
    data = sys.argv[1]
    cg_props = sys.argv[2]
    name = sys.argv[3]
    compound_props = json.load(open(cg_props))
    ml_convert = json.load(open("mlConvert.json"))
    with open(data) as data_source, open(name + ".csv", "w") as data_out:
        reader = csv.reader(data_source)
        writer = csv.writer(data_out)
        reader.next()
        writer.writerow(rxn_calculator.headers)
        rows = [line for line in reader]
        #prior = build_prior(rows)
        reader = rows# + prior
        for line in reader:
            try:
                writer.writerow(parse_rxn(line, compound_props, ml_convert))
            except Exception as e:
                print "exception", e, line


if __name__ == "__main__":
    main()
