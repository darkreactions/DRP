import sys, subprocess
import rdkit.Chem as Chem

class CGCalculator:
    def __init__(self, compound, file_name, smiles, m_type,jchem_path = "/home/praccugl/ChemAxon/JChem/bin", sdf_path="sdf", error_pipe = sys.stderr):
        ''' smiles = SMILES string of the compound
        m_type = the type of the molecule. Can be: 
            'w' (water),'p' (pH),'i' (inorganic),'o' (organic),'ox' (oxalate)
        '''

        self.verbose = False

        self.sdf_path = sdf_path
        self.jchem_path = jchem_path
        self.cxcalc = "{0}/cxcalc -N ih".format(self.jchem_path)
        self.error_pipe = error_pipe

        self.m_type = m_type
        self.smiles = smiles 
        self.compound = compound
        self.file_name = file_name

        self.sdf = "{0}/{1}".format(self.sdf_path, file_name) 
        self.calc_sdf()

        self.properties_map = {"mw": self.calc_MW(), "NopH": self.calc_NopH(), "polsurf": self.calc_polSurf(), "msacc": self.calc_msacc(), "msdon": self.calc_msdon(), "projectionArea": self.calc_projectionArea(), "type": self.m_type, "atoms": list(set(atoms_from_smiles(smiles))), "compound":compound, "smiles": self.smiles}

    def get_properties(self):
        return self.properties_map

    def calc_sdf(self):
        ''' JChem's molconvert takes the smiles string and generates
        an SDF-format string. We can use this string to calculate all the
        properties we are interested in. This method currently doesn't error
        check'''
        run = "{0}/molconvert \"sdf\" -3:S{{nofaulty}} -s ".format(self.jchem_path) 
        run += "\"{0}\" > {1}".format(self.smiles,self.sdf)
        if self.verbose: print run
        return subprocess.check_output(run, shell=True)


    def calc_MW(self):
        ''' Calculates the molecular weight of the compound'''
        run = "{0} mass \"{1}\"".format(self.cxcalc, self.sdf)
        if self.verbose: print run
        mw = subprocess.check_output(run, shell=True)
        #if self.abbrev == "V2O5": # hack to double the moles for V2
        #    mw = float(mw)/2.0 # TODO: figure out why we did this and how to generalize
        return float(mw)

    def calc_NopH(self):
        ''' Calculates 'pH insensitive descriptors' '''
        if False: #self.m_type == "Inorg" or self.m_type == "pH":
            return [-1,-1]
        run = "{0} avgpol refractivity {1}".format(self.cxcalc, self.sdf)
        if self.verbose: print run
        nopH = subprocess.check_output(run, shell=True)

        nopH_properties = [float(prop) for prop in nopH.split()]
        return nopH_properties

    def calc_polSurf(self):
        ''' This shouldn't be called polSurf, but it is...
        This actually calculates the pH sensitive properties for each pH
        The output is a list of length 9*14. For every pH in [1 ..14],
        we are calculating 9 properties. So the first 9 values are the
        properties for pH 1. The next 9 are for pH 2. etc.
        '''

        if self.m_type == "Inorg" or self.m_type == "pH":
            return [-1 for i in range(126)]
        run = self.cxcalc + " "
        for pH in range(1,15):
            run += "avgpol -H {0} molpol -H {0} msa -t ".format(pH)
            run += "\"vanderwaals,ASA,ASA+,ASA-,ASA_H,ASA_P\" "
            run += "-H {0} polarsurfacearea -H {0} ".format(pH)
        run += " {0}".format(self.sdf)
        if self.verbose: print run
        raw = subprocess.check_output(run, shell=True)
        pH_properties = [float(prop) for prop in raw.rstrip().split()]
        return pH_properties
        

    def calc_msacc(self):
        if self.m_type == "Inorg" or self.m_type == "pH":
            return [-1 for i in range(15)]
        ''' TODO: some kind of acceptor '''
        run = "{0} hbda -t msacc {1}".format(self.cxcalc, self.sdf)
        if self.verbose: print run
        raw = subprocess.check_output(run, shell=True)
        msacc = [float(prop) for prop in raw.rstrip().split()]
        return msacc

    def calc_msdon(self):
        ''' TODO: some kind of donor '''
        if self.m_type == "Inorg" or self.m_type == "pH":
            return [-1 for i in range(15)]

        run = "{0} hbda -t msdon {1}".format(self.cxcalc, self.sdf)
        if self.verbose: print run
        raw = subprocess.check_output(run, shell=True)
        msdon = [float(prop) for prop in raw.rstrip().split()] 
        return msdon

    def make_nonsalt_sdf(self):
        fragments = self.smiles.split(".")
        fragments = [frag for frag in fragments if not is_salt(frag)]
        if self.verbose: print fragments

        if len(fragments) > 1:
            raise Exception("calc_projectionArea > make_nonsalt_sdf > Failed to remove excess fragments from {0}".format(self.smiles))
        frag = fragments[0] 
        run = "{0}/molconvert \"sdf\" -3:S{{nofaulty}} -s ".format(self.jchem_path)
        run += "\"{0}\" > {1}.nonsalt".format(frag, self.sdf )
        if self.verbose: print run
    	subprocess.check_output(run, shell=True)
        return self.sdf + ".nonsalt" 

    def calc_projectionArea(self):
        ''' TODO: calculate projection area descriptors. Note that the
        fragmentcount has to be 1 for this to be defined. This is not
        the case for some of the inorganics or for dihydrochlorides. 
        For the latter, we can do a substitution to get this. Note that 
        these are also pH insensitive.'''

        failure = [-1 for i in range(6)] 
        if self.m_type == "Inorg" or self.m_type == "pH":
            return failure
        sdf = self.sdf
        run = "{0} fragmentcount {1}".format(self.cxcalc, sdf)
        if self.verbose: print run
        num_frag = subprocess.check_output(run, shell=True).rstrip()

        if int(num_frag) != 1: 
            if "." in self.smiles: # there are fragments
                try:
                    sdf = self.make_nonsalt_sdf()
                except Exception as e:
                    #self.error_pipe.write(str(e))
                    print str(e)
                    return failure
            else:
                sdf = self.sdf
    
            run = "{0} fragmentcount {1}".format(self.cxcalc, sdf)
            if self.verbose: print run
            num_frag = subprocess.check_output(run, shell=True).rstrip()
            
            if int(num_frag) != 1:
                self.error_pipe.write("\ncalc_projectionArea > Didn't find lone fragment in {0}\n".format(self.smiles))
                return failure

        run = "{0} maximalprojectionarea maximalprojectionradius".format(self.cxcalc)
        run += " maximalprojectionsize minimalprojectionarea"
        run += " minimalprojectionradius minimalprojectionsize"
        run += " {0}".format(sdf)

        if self.verbose: print run
        raw_out = subprocess.check_output(run, shell=True).rstrip()
        return [float(prop) for prop in raw_out.split()]
        

def is_salt(fragment_smiles):
    def match_sulfate(atom_list):
        return len(atom_list) == 5 and atom_list.count('o') == 4 and atom_list.count('s') == 1

    def match_ammonium_meta(atom_list):
        return (len(atom_list) == 5 and atom_list.count('h') == 4 and atom_list.count('n') == 1) or (len(atom_list) == 1 and atom_list.count('n') == 1)

    def match_chloride(atom_list):
        return (len(atom_list) == 1 and atom_list.count('cl') == 1)

    def match_iodine(atom_list):
        return (len(atom_list) == 1 and atom_list.count('i') == 1)

    def match_bromide(atom_list):
        return (len(atom_list) == 1 and atom_list.count('br') == 1)
    def match_fluorine(atom_list):
        return (len(atom_list)) == 1 and atom_list.count('f') == 1

    matcher_list = [ match_sulfate, match_ammonium_meta, match_chloride, match_iodine, match_bromide, match_fluorine ]

    atom_list = [atom.lower() for atom in atoms_from_smiles(fragment_smiles)]

    for matcher in matcher_list:
        if matcher(atom_list):
            return True
    return False


def atoms_from_smiles(smiles):
    mols = Chem.MolFromSmiles(smiles,sanitize=False)
    if mols == None:
        return []
    atoms = mols.GetAtoms()
    return [atom.GetSymbol() for atom in atoms]
