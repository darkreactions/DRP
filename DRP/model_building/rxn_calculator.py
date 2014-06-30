import sys, os, subprocess
import json, csv, math
import rdkit.Chem as Chem

django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.settings import RESEARCH_DIR

class PropertiesCalculator:
    def __init__(self, compoundMoles, inorganicList, compoundAcc, compoundDon, organicList, isWater):
        self.compoundMoles = compoundMoles
        self.inorganicList = inorganicList
        self.organicList = organicList
        self.compoundAcc = compoundAcc
        self.compoundDon = compoundDon
        self.isWater = isWater

    #inorganic/water mole ratio
    def inorg_water_mole_ratio(self):
        totalInorganic = 0.0
        for k in range(5):
            totalInorganic += self.compoundMoles[k]*self.inorganicList[k]
        totalWater = self.compoundMoles[self.isWater]
        #if totalWater == 0:
        #    totalWater = 0.00001
        return totalInorganic/totalWater


    # organic / water mole ratio
    def org_water_mole_ratio(self):
        totalOrganic = 0.0
        for k in range(5):
            totalOrganic += self.compoundMoles[k]*self.organicList[k]
        totalWater = self.compoundMoles[self.isWater]
        #if totalWater == 0:
        #    totalWater = 0.00001

        return totalOrganic/totalWater

    def org_water_acceptor_on_inorg_ratio(self):
        #organic/water acceptor-on-inorganic ratio
        #TODO: the name doesn't make sense
        totalOrganic = 0.0
        for k in range(5):
            totalOrganic += self.compoundMoles[k]*self.organicList[k]*self.compoundAcc[k]
        return totalOrganic/(self.compoundMoles[self.isWater]*self.compoundDon[self.isWater])

    def org_water_donor_on_inorg_ratio(self):
        #organic/water donor-on-inorganic ratio
        # the name doesn't make sense
        totalOrganic = 0.0
        for k in range(5):
            totalOrganic += self.compoundMoles[k]*self.organicList[k]*self.compoundDon[k]
        return totalOrganic/(self.compoundMoles[self.isWater]*self.compoundAcc[self.isWater])

    def inorg_org_mole_ratio(self):
        #inorganic/organic mole-ratio
        totalInorganic = 0.0
        totalOrganic = 0.0
        for k in range(5):
            totalOrganic += self.compoundMoles[k]*self.organicList[k]
            totalInorganic += self.compoundMoles[k]*self.inorganicList[k]
        if not totalOrganic:
            totalOrganic = 0.00001
        return totalInorganic/totalOrganic

    def not_water_mole_ratio(self):
        #not water/water mole ratio
        totalInorganic = 0.0
        totalOrganic = 0.0
        for k in range(5):
            totalOrganic += self.compoundMoles[k]*self.organicList[k]
            totalInorganic += self.compoundMoles[k]*self.inorganicList[k]
        totalWater = self.compoundMoles[self.isWater]
        if totalWater == 0:
            totalWater = 0.00001
        return (totalInorganic+totalOrganic)/totalWater


def fixmL(t):
    o = ''
    for i in t:
        if ord(i) == ord('/') or ord(i) > ord('9') or ord(i) < ord('.'):
            return o
        o += i
    return o
def dictFix(d):
    return dict(((key, float(d[key])) for key in d))

def distList(indicator, properties_lists):
    from operator import mul
    ''' Every item in "properties_lists" is a list of values for a set of
    properties. Alternatively: properties_list is an n x 19 matrix, with 19
    features, where each row corrresponds to a compound.

    Only some rows are relevant, so we filter out those rows that indicator
     does not carry a 1 for....
    For each column, I want to take the min, max, average and geometric average
    '''

    def min_f(l):
        if len(l) == 0: return 0
        return min(l)

    def max_f(l):
        if len(l) == 0: return 0
        return max(l)

    def mean(l):
        if len(l) == 0: return 0
        return sum(l)/float(len(l))

    def geoavg(l):
        if len(l) == 0: return 0
        return reduce(mul, l)**(1.0/float(len(l)))

    assert(len(properties_lists) == len(indicator))
    compressed = [properties_lists[i] for i in range(len(indicator)) if indicator[i] == 1]
    transposed = zip(*compressed) # each element is a list of values for that column's property
    filtered = map(lambda property_set: filter(lambda prop_val: prop_val != -1, property_set), transposed)

    minned = map(min_f, filtered)
    maxxed = map(max_f, filtered)
    meaned = map(mean, filtered)
    geomed = map(geoavg, filtered)
    if len(minned) == 0:
        raise Exception("empty distlist")
    return maxxed +  minned + meaned + geomed


def atoms_from_smiles(smiles):
    mols = Chem.MolFromSmiles(str(smiles),sanitize=False)
    if mols == None:
        return []
    atoms = mols.GetAtoms()
    return [atom.GetSymbol() for atom in atoms]

def atomic_properties(atom_list, smiles_pairs, counts = None):
    interesting = ["Se", "Te", "V", "Mo", "Zn","Ga","Co", "Cr"]
    if not counts:
        atom_count = [(filter(lambda x: x in interesting, atoms_from_smiles(k[0])), k[1]) for k in smiles_pairs]

        counts = {}
        for p in atom_count:
            for a in p[0]:
                if a not in counts:
                    counts[a] = 0
                counts[a] += min(1,p[0].count(a))*p[1]


    atoms = {atom: properties[atom] for atom in atom_list if atom in interesting}

    props = []

    props.append("yes" if any([atoms[b]["Actinide"] for b in atoms]) else "no")
    props.append("yes" if any([atoms[b]["AlkaliMetal"] for b in atoms]) else "no")
    props.append("yes" if any([atoms[b]["Lanthanide"] for b in atoms]) else "no")

    for i in range(1,8):
        props.append("yes" if any([atoms[b]["Period"] == i for b in atoms]) else "no")

    for i in range(1,19):
        props.append("yes" if any([atoms[b]["Group"] == i for b in atoms]) else "no")

    for i in range(0,8):
        props.append("yes" if any([atoms[b]["Valence"] == i for b in atoms]) else "no")

    numeric_props = [ [ atoms[b]["IonizationEnergies"], atoms[b]["ElectronAffinity"], atoms[b]["PaulingElectronegativity"], atoms[b]["PearsonElectronegativity"], atoms[b]["hardness"], atoms[b]["AtomicRadius"] ] for b in atoms]
    numeric_props_weighted = [ map(lambda x: x*counts[b],[ atoms[b]["IonizationEnergies"], atoms[b]["ElectronAffinity"], atoms[b]["PaulingElectronegativity"], atoms[b]["PearsonElectronegativity"], atoms[b]["hardness"], atoms[b]["AtomicRadius"] ]) for b in atoms]
    indicator = [1 for i in range(len(atoms))]
    try:
        props += distList(indicator, numeric_props)
        props += distList(indicator, numeric_props_weighted)
    except Exception as e:
	err = str(atom_list) + " not interesting?"
        #raise Exception()
	print err
	props += [0 for i in range(48)]
    return props

properties = json.load(open(RESEARCH_DIR + "scripts/atomic_props.json"))


bools = ["Actinide", "AlkaliMetal", "Lanthanide",
    "P1", "P2", "P3", "P4", "P5", "P6", "P7",
    "G1", "G2", "G3", "G4", 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16', 'G17', 'G18',
    "V0", "V1", "V2", "V3", "V4", "V5", "V6", "V7"]

field_names = ["Actinide", "AlkaliMetal", "Lanthanide",
    "P1", "P2", "P3", "P4", "P5", "P6", "P7",
    "G1", "G2", "G3", "G4", 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16', 'G17', 'G18',
    "V0", "V1", "V2", "V3", "V4", "V5", "V6", "V7",
    "IonizationMax", "EAMax", "PaulingElectronegMax", "PearsonElectronegMax", "hardnessMax", "AtomicRadiusMax",
    "IonizationMin", "EAMin", "PaulingElectronegMin", "PearsonElectronegMin", "hardnessMin", "AtomicRadiusMin",
    "IonizationMean", "EAMean", "PaulingElectronegMean", "PearsonElectronegMean", "hardnessMean", "AtomicRadiusMean",
    "IonizationGeom", "EAGeom", "PaulingElectronegGeom", "PearsonElectronegGeom", "hardnessGeom", "AtomicRadiusGeom",
    "IonizationMaxWeighted", "EAMaxWeighted", "PaulingElectronegMaxWeighted", "PearsonElectronegMaxWeighted", "hardnessMaxWeighted", "AtomicRadiusMaxWeighted",
    "IonizationMinWeighted", "EAMinWeighted", "PaulingElectronegMinWeighted", "PearsonElectronegMinWeighted", "hardnessMinWeighted", "AtomicRadiusMinWeighted",
    "IonizationMeanWeighted", "EAMeanWeighted", "PaulingElectronegMeanWeighted", "PearsonElectronegMeanWeighted", "hardnessMeanWeighted", "AtomicRadiusMeanWeighted",
    "IonizationGeomWeighted", "EAGeomWeighted", "PaulingElectronegGeomWeighted", "PearsonElectronegGeomWeighted", "hardnessGeomWeighted", "AtomicRadiusGeomWeighted",
    ]
atomsz = ['Na', 'Li', 'Te', 'Br', 'K', 'C', 'F', 'I', 'Mo', 'O', 'N', 'P', 'S', 'V', 'Se', 'Zn', 'Co', 'Cl', 'Ga', 'Cs', 'Cr', 'Cu']
headers = ['XXXtitle', 'XXXinorg1', 'XXXinorg1mass',
            'XXXinorg1moles', 'XXXinorg2', 'XXXinorg2mass',
            'XXXinorg2moles', 'XXXinorg3', 'XXXinorg3mass','XXXinorg3moles', 'XXXorg1', 'XXXorg1mass',
            'XXXorg1moles', 'XXXorg2', 'XXXorg2mass',
            'XXXorg2moles', 'XXXoxlike1', 'XXXoxlike1mass',
            'XXXoxlike1moles', 'Temp_max', 'time', 'slowCool', 'pH',
            'leak', 'numberInorg', 'numberOrg', 'numberOxlike',
            'numberComponents', 'orgavgpolMax', 'orgrefractivityMax',
            'orgmaximalprojectionareaMax',
            'orgmaximalprojectionradiusMax',
            'orgmaximalprojectionsizeMax',
            'orgminimalprojectionareaMax',
            'orgminimalprojectionradiusMax',
            'orgminimalprojectionsizeMax',
            'orgavgpol_pHdependentMax', 'orgmolpolMax',
            'orgvanderwaalsMax', 'orgASAMax', 'orgASA+Max',
            'orgASA-Max', 'orgASA_HMax', 'orgASA_PMax',
            'orgpolarsurfaceareaMax', 'orghbdamsaccMax',
            'orghbdamsdonMax', 'orgavgpolMin', 'orgrefractivityMin',
            'orgmaximalprojectionareaMin',
            'orgmaximalprojectionradiusMin',
            'orgmaximalprojectionsizeMin',
            'orgminimalprojectionareaMin',
            'orgminimalprojectionradiusMin',
            'orgminimalprojectionsizeMin', 'orgavgpol_pHdependentMin',
            'orgmolpolMin', 'orgvanderwaalsMin', 'orgASAMin',
            'orgASA+Min', 'orgASA-Min', 'orgASA_HMin', 'orgASA_PMin',
            'orgpolarsurfaceareaMin', 'orghbdamsaccMin',
            'orghbdamsdonMin', 'orgavgpolArithAvg',
            'orgrefractivityArithAvg',
            'orgmaximalprojectionareaArithAvg',
            'orgmaximalprojectionradiusArithAvg',
            'orgmaximalprojectionsizeArithAvg',
            'orgminimalprojectionareaArithAvg',
            'orgminimalprojectionradiusArithAvg',
            'orgminimalprojectionsizeArithAvg',
            'orgavgpol_pHdependentArithAvg',
            'orgmolpolArithAvg', 'orgvanderwaalsArithAvg',
            'orgASAArithAvg', 'orgASA+ArithAvg', 'orgASA-ArithAvg',
            'orgASA_HArithAvg', 'orgASA_PArithAvg',
            'orgpolarsurfaceareaArithAvg', 'orghbdamsaccArithAvg',
            'orghbdamsdonArithAvg', 'orgavgpolGeomAvg',
            'orgrefractivityGeomAvg',
            'orgmaximalprojectionareaGeomAvg',
            'orgmaximalprojectionradiusGeomAvg',
            'orgmaximalprojectionsizeGeomAvg',
            'orgminimalprojectionareaGeomAvg',
            'orgminimalprojectionradiusGeomAvg',
            'orgminimalprojectionsizeGeomAvg',
            'orgavgpol_pHdependentGeomAvg',
            'orgmolpolGeomAvg', 'orgvanderwaalsGeomAvg',
            'orgASAGeomAvg', 'orgASA+GeomAvg', 'orgASA-GeomAvg',
            'orgASA_HGeomAvg', 'orgASA_PGeomAvg',
            'orgpolarsurfaceareaGeomAvg', 'orghbdamsaccGeomAvg',
            'orghbdamsdonGeomAvg', 'oxlikeavgpolMax',
            'oxlikerefractivityMax', 'oxlikemaximalprojectionareaMax',
            'oxlikemaximalprojectionradiusMax',
            'oxlikemaximalprojectionsizeMax',
            'oxlikeminimalprojectionareaMax',
            'oxlikeminimalprojectionradiusMax',
            'oxlikeminimalprojectionsizeMax'] + [
            'oxlikeavgpol_pHdependentMax', 'oxlikemolpolMax',
            'oxlikevanderwaalsMax', 'oxlikeASAMax', 'oxlikeASA+Max',
            'oxlikeASA-Max', 'oxlikeASA_HMax', 'oxlikeASA_PMax',
            'oxlikepolarsurfaceareaMax', 'oxlikehbdamsaccMax',
            'oxlikehbdamsdonMax', 'oxlikeavgpolMin',
            'oxlikerefractivityMin', 'oxlikemaximalprojectionareaMin',
            'oxlikemaximalprojectionradiusMin',
            'oxlikemaximalprojectionsizeMin',
            'oxlikeminimalprojectionareaMin',
            'oxlikeminimalprojectionradiusMin',
            'oxlikeminimalprojectionsizeMin',
            'oxlikeavgpol_pHdependentMin', 'oxlikemolpolMin',
            'oxlikevanderwaalsMin', 'oxlikeASAMin', 'oxlikeASA+Min',
            'oxlikeASA-Min', 'oxlikeASA_HMin', 'oxlikeASA_PMin',
            'oxlikepolarsurfaceareaMin', 'oxlikehbdamsaccMin',
            'oxlikehbdamsdonMin', 'oxlikeavgpolArithAvg',
            'oxlikerefractivityArithAvg',
            'oxlikemaximalprojectionareaArithAvg',
            'oxlikemaximalprojectionradiusArithAvg',
            'oxlikemaximalprojectionsizeArithAvg',
            'oxlikeminimalprojectionareaArithAvg',
            'oxlikeminimalprojectionradiusArithAvg',
            'oxlikeminimalprojectionsizeArithAvg',
            'oxlikeavgpol_pHdependentArithAvg',
            'oxlikemolpolArithAvg', 'oxlikevanderwaalsArithAvg',
            'oxlikeASAArithAvg', 'oxlikeASA+ArithAvg',
            'oxlikeASA-ArithAvg', 'oxlikeASA_HArithAvg',
            'oxlikeASA_PArithAvg', 'oxlikepolarsurfaceareaArithAvg',
            'oxlikehbdamsaccArithAvg', 'oxlikehbdamsdonArithAvg',
            'oxlikeavgpolGeomAvg', 'oxlikerefractivityGeomAvg',
            'oxlikemaximalprojectionareaGeomAvg',
            'oxlikemaximalprojectionradiusGeomAvg',
            'oxlikemaximalprojectionsizeGeomAvg',
            'oxlikeminimalprojectionareaGeomAvg',
            'oxlikeminimalprojectionradiusGeomAvg',
            'oxlikeminimalprojectionsizeGeomAvg',
            'oxlikeavgpol_pHdependentGeomAvg', 'oxlikemolpolGeomAvg',
            'oxlikevanderwaalsGeomAvg', 'oxlikeASAGeomAvg',
            'oxlikeASA+GeomAvg', 'oxlikeASA-GeomAvg',
            'oxlikeASA_HGeomAvg', 'oxlikeASA_PGeomAvg',
            'oxlikepolarsurfaceareaGeomAvg', 'oxlikehbdamsaccGeomAvg',
            'oxlikehbdamsdonGeomAvg', 'inorg-water-moleratio',
            'org-water-moleratio', 'orgacc-waterdonratio',
            'orgdon-wateraccratio', 'inorg-org-moleratio',
            'notwater-water-moleratio'] + atomsz + field_names + ['purity', 'outcome']
