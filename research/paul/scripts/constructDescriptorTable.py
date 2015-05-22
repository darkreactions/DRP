#!/usr/bin/python
# -*- coding: latin-1 -*-
import sys, os, subprocess
import json, csv, math
import traceback

def next_alternate(abbrev):
    if abbrev[-1] == 'z':
        return abbrev + 'a'
    else:
        return abbrev[:-1] + chr(ord(abbrev[-1])+1)

def fixmL(t):
    o = ''
    for i in t:
        if ord(i) == ord('/') or ord(i) > ord('9') or ord(i) < ord('.'):
            return o
        o += i
    return o
def dictFix(d):
    return dict(((key, float(d[key])) for key in d))

def distList(categoryList, propertyList):
    minList, maxList,geomList,arithList,foo = ([],[],[],[],[])
    for i in range(5):
        if categoryList[i] == 1:
            foo = propertyList[i]

            if len(minList):
                for j in range(len(foo)):
                    if foo[j] < 0:
                        continue
                    minList[j] = min(minList[j], foo[j])
                    maxList[j] = max(maxList[j], foo[j])
                    geomList[j] *= foo[j]
                    arithList[j] += foo[j]
            else:
                minList = foo[:]
                maxList = foo[:]
                if -1 in foo[:]:
                    geomList = [abs(x) for x in foo] 
                    arithList = [0 for thing in foo]
                else:
                    geomList = foo[:]
                    arithList = foo[:]
    return maxList + minList + arithList + geomList

def fault_args(errstr):
    raise Exception('''Error: %s \n\nUsage:
    contructDescriptorTable.py 
    <compound_guide_path>
    <reaction_list_path> <out_file> 
    [--testCompoundGuide|--testFragment|--existingCG|--validateOnly]''' % errstr)


if __name__ == "__main__":
    v_tot = 0.023
    p_init = 1.00
    t_init = 273.0 + 20.0
    R = 0.0831451
    rho_water = 0.001

    z_b_shift = 14 # how many places to shift for dbz
    # HANDLE ARGUMENTS
    if len(sys.argv) < 4:
        fault_args('Insufficient arguments!')
        
    
    test_cg = False
    test_frag = False
    test_existing = False 
    validate_only = False
    
    if len(sys.argv) > 4:
        try:
            #Clean and strip the optional argument.
            option_arg = str(sys.argv[4]).lower().replace("-","").replace("_","")
            
            if option_arg == "testcompoundguide":
                test_cg = True 
            elif option_arg == "testfragment":
                test_frag = True 
            elif option_arg == "existingcg":
                test_existing = True
            elif option_arg == "validateonly":
                validate_only = True
            elif len(sys.argv) > 5:
                raise Exception('More than 4 arguments?')
            else:
                raise Exception('Unknown fourth argument!')            
        except Exception as e:
            fault_args(e)

    ########## Set up file paths ######################################
    #Translate arguments.
    rxn_database_path = sys.argv[1]
    compound_guide_path = sys.argv[2]
    prefix = sys.argv[3]
    
    #Assume that script is in scripts folder.
    CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    TRIAL_DIR = "{0}/{1}".format(CHEMML_DIR, prefix)
    SDF_DIR = "{0}/sdf".format(TRIAL_DIR)
    
    #Create new directories if they do not already exist.
    if not os.path.exists(TRIAL_DIR):
        os.makedirs(TRIAL_DIR)
    if not os.path.exists(SDF_DIR):
        os.makedirs(SDF_DIR)
    
    # Get JChem directory from config.json:
    with open(CHEMML_DIR+"/config/config.json", "r") as CONFIG_FILE:
        jchem_path = json.load(CONFIG_FILE)["JCHEM_DIRECTORY"]
    
    
    mlConvert_json_path = "{0}/scripts/mlConvert.json".format(CHEMML_DIR)
    restart_json_path = "{0}/restart.json".format(TRIAL_DIR, prefix)
    smiles_json_path = "{0}/smiles.json".format(TRIAL_DIR, prefix)
    
    #Create the hypothetical CSV file instead if the actual expanded CSV exists.
    expanded_csv_path = "{0}/{1}_expanded.csv".format(TRIAL_DIR, prefix)
    
    abbrev_to_IUPAC = {}
    abbrev_to_type = {}
    abbrev_to_MW = {}
    abbrev_to_NopH = {}
    abbrev_to_polSurf = {}
    abbrev_to_msacc = {}
    abbrev_to_msdon = {}
    abbrev_to_nonsaltIUPAC = {}
    abbrev_to_ProjectionArea = {}
    IUPAC_to_SMILES = {}
    skip = {}
    
    if test_existing:
        existing_dict = json.load(open(restart_json_path))
        abbrev_to_IUPAC = existing_dict['abbrev_to_IUPAC']
        abbrev_to_type = existing_dict['abbrev_to_type']
        abbrev_to_MW = existing_dict['abbrev_to_MW']
        abbrev_to_NopH = existing_dict['abbrev_to_NopH']
        abbrev_to_polSurf = existing_dict['abbrev_to_polSurf']
        abbrev_to_msacc = existing_dict['abbrev_to_msacc']
        abbrev_to_msdon = existing_dict['abbrev_to_msdon']
        abbrev_to_nonsaltIUPAC = existing_dict['abbrev_to_nonsaltIUPAC']
        abbrev_to_ProjectionArea = existing_dict['abbrev_to_ProjectionArea']
        skip = existing_dict['skip']
    else:
        with open(compound_guide_path) as guide_file:
            alternate_abbrev = 'a'
            guide_file.next() # discard comment line? Per JS's code
            for row in csv.reader(guide_file):
                abbrev = ''.join(row[0].split()) # Remove whitespace
                IUPAC = row[1].replace('’', "'")
                try:
                    if not abbrev or not IUPAC: ###NEEDED?
                        sys.stderr.write("WARNING: no abbrev or IUPAC! {} , {}? SKIPPING.\n".format(abbrev, IUPAC)) 
                        continue
                    m_type = row[3].lower() # named type in Josh Schrier's script 
                    # this is named "m_type" because "type" shadows the python "type" 
                    m_type = m_type.replace('inorg','i').replace('org','o').replace('ph','p').replace('oxalate','ox').replace('water','w').replace('sol','s')

                    if m_type in ['s','p']:
                        skip[abbrev] = m_type
                        continue
                    if m_type not in ['w','i','o','ox']:
                        err = "WARNING: invalid reactant type, --%s" % m_type
                        err += "--specified for species %s = " % abbrev
                        err += "%s .. SKIPPING\n" % IUPAC
                        sys.stderr.write(err)
                        continue
                    abbrev_to_type[abbrev] = m_type
                    abbrev_to_IUPAC[abbrev] = IUPAC


                    try:
                        run = jchem_path + "/molconvert smiles:n -s \"%s\"" % IUPAC
                        smiles = subprocess.check_output(run, shell=True).rstrip() 
                        IUPAC_to_SMILES[IUPAC] = smiles.split()[0]
                    except Exception as e:
                        sys.stderr.write("WARNING: failed to make smiles for %s: %s\n" % (IUPAC, str(e)))

                    if validate_only:
                        sys.stderr.write(IUPAC + "=" + abbrev + "\n")
                        continue # don't generate fragments for validation

                    ################# Generate SDF files ##############

                    abbrev_filename = abbrev
                    if not abbrev_filename.isalnum():
                        abbrev_filename = alternate_abbrev
                        alternate_abbrev = next_alternate(alternate_abbrev)
                    run = "{0}/molconvert \"sdf\" -3:S{{nofaulty}} -s ".format(jchem_path)
                    run += "\"{0}\" > {1}/{2}".format(IUPAC, SDF_DIR, abbrev_filename)
                    subprocess.call(run, shell=True)

                    # note from JS: below is the most readable. If we need to
                    # improve performance, it is more efficient to combine
                    # all of the cxcalc executions into one and then split
                    # up the resulting string within perl.
                    if not test_frag:
                        # compute molecular weight data
                        run = "{0}/cxcalc -N ih mass {1}/{2}".format(jchem_path,SDF_DIR, abbrev_filename)
                        sys.stderr.write(run+"\n")
                        mw = subprocess.check_output(run, shell=True).rstrip()
                        if abbrev == "V2O5": #hack to double the moles for V2
                            mw = float(mw)/2.0
                        abbrev_to_MW[abbrev] = float(mw)
                        if test_cg:
                            sys.stderr.write("{} --> {}\n".format(abbrev, mw))
                        # compute pH insensitive descriptors
                        run = "{0}/cxcalc -N ih avgpol refractivity {1}/{2}".format(jchem_path, SDF_DIR, abbrev_filename)
                        nopH = subprocess.check_output(run, shell=True).rstrip()
                        
                        abbrev_to_NopH[abbrev] = [float(o) for o in nopH.split()]
                        if test_cg:
                            sys.stderr.write("{} --> {}\n".format(abbrev,abbrev_to_NopH[abbrev])) 

                        # compute pH sensitive descriptor sets
                        # except for # hbond donor acceptors
                        run = "%s/cxcalc -N ih " % jchem_path 
                        for pH in range(1,15):
                            run += "avgpol -H %d molpol -H %d msa -t " % (pH, pH) 
                            run += "\"vanderwaals,ASA,ASA+,ASA-,ASA_H,ASA_P\" "
                            run += "-H %d polarsurfacearea -H %d " % (pH, pH)
                        run += " {0}/{1}".format(SDF_DIR, abbrev_filename) 
                        abbrev_to_polSurf[abbrev] = map(float, subprocess.check_output(
                            run, shell=True).rstrip().split()) #TODO: fix jank
                        #elements (pH-1)*9 .. (pH*9-1) in this string are for the
                        # reference pH (i.e., 9 elements per pH value)

                        if test_cg:
                            sys.stderr.write("%s --> %s\n" % (abbrev, 
                                    abbrev_to_polSurf[abbrev]))

                        run = "{0}/cxcalc -N ih hbda -t msacc {1}/{2}".format(jchem_path, SDF_DIR, abbrev_filename) 
                        abbrev_to_msacc[abbrev] = map(float, subprocess.check_output(run,
                                shell=True).rstrip().split())

                        run = "{0}/cxcalc -N ih hbda -t msdon {1}/{2}".format(jchem_path, SDF_DIR, abbrev_filename) 
                        abbrev_to_msdon[abbrev] = map(float, subprocess.check_output(
                                run, shell=True).rstrip().split())
                        if test_cg:
                            sys.stderr.write("{} --> {}\n".format(abbrev, abbrev_to_msdon[abbrev]))

                    # compute projection area descriptors. Note that the 
                    # fragmentcount has to be 1 for this to be defined. This
                    # is not the case for some of the inorganics or for the
                    # dihydrochlorides. For the latter, we can do a
                    # substition to get this. Note that these are also pH
                    # insensitive.
                    n_frag = 0
                    nonsalt_IUPAC = ''
                    if 'dihydrochloride' in IUPAC:
                        nonsalt_IUPAC = IUPAC.replace('dihydrochloride','')
                    elif 'ammonium meta' in IUPAC:
                        nonsalt_IUPAC = IUPAC.replace('ammonium meta','')
                    elif ' sulfate' in IUPAC:
                        nonsalt_IUPAC = IUPAC.replace(' sulfate','')
                    if nonsalt_IUPAC:
                        # create a substitute version to compute the
                        # project area descriptors
                        abbrev_to_nonsaltIUPAC[abbrev] = nonsalt_IUPAC
                        # compute the fragment number for the non-salt
                        if True: #not os.path.exists("sdf/%s.nonsalt" % abbrev):
                            run = "%s/molconvert \"sdf\" -3:S{nofaulty} -s " % jchem_path 
                            run += "\"{0}\" > {1}/{2}.nonsalt".format(nonsalt_IUPAC, SDF_DIR, abbrev_filename)
                            subprocess.call(run, shell=True)
                            run = "{}/cxcalc -N ih fragmentcount {}/{}.nonsalt".format(jchem_path, SDF_DIR, abbrev_filename)
                    else:
                        if True: # not os.path.exists("sdf/%s" % abbrev):
                            run = "{}/cxcalc -N ih fragmentcount {}/{}".format(jchem_path, SDF_DIR, abbrev_filename)
                    n_frag = subprocess.check_output(run, shell=True).rstrip()
                    if int(n_frag) == 1: # compute projection area descirptors
                        run = "%s/cxcalc -N ih maximalprojectionarea maximalpro" % jchem_path
                        run += "jectionradius maximalprojectionsize minimalpro"
                        run += "jectionarea minimalprojectionradius minimalpro"
                        run += "jectionsize {}/{}".format(SDF_DIR, abbrev_filename)
                        if abbrev in abbrev_to_nonsaltIUPAC.keys():
                            run += ".nonsalt"
                        foo = subprocess.check_output(run, shell=True).rstrip()
                        if test_frag:
                            if abbrev_to_nonsaltIUPAC[abbrev]:
                                sys.stderr.write("Warning: nonsalt.\n")
                            sys.stderr.write("%s --> %s\n" % (abbrev, foo))
                        abbrev_to_ProjectionArea[abbrev] = map(float, foo.split())
                    else:
                        run = "%s = %s has %s fragments" % (abbrev, IUPAC, n_frag)
                        run += "...cannot compute projection area"
                        abbrev_to_ProjectionArea[abbrev] = [-1,-1,-1,-1,-1,-1]
                except Exception as e:
                    sys.stderr.write("Merp? EXCEPTION: %s\n" % str(e))
                    if abbrev in abbrev_to_IUPAC.keys():
                        del abbrev_to_IUPAC[abbrev]
                    if abbrev in abbrev_to_type.keys():
                        del abbrev_to_type[abbrev]
                    if abbrev in abbrev_to_MW.keys():
                        del abbrev_to_MW[abbrev]
                    if abbrev in abbrev_to_NopH.keys():
                        del abbrev_to_NopH[abbrev]
                    if abbrev in abbrev_to_polSurf.keys():
                        del abbrev_to_polSurf[abbrev]
                    if abbrev in abbrev_to_msacc.keys():
                        del abbrev_to_msacc[abbrev]
                    if abbrev in abbrev_to_msdon.keys():
                        del abbrev_to_msdon[abbrev]
                    if abbrev in abbrev_to_nonsaltIUPAC.keys():
                        del abbrev_to_nonsaltIUPAC[abbrev]
                    if abbrev in abbrev_to_ProjectionArea.keys():
                        del abbrev_to_ProjectionArea[abbrev]
                    if IUPAC in IUPAC_to_SMILES.keys():
                        del IUPAC_to_SMILES[IUPAC]
                    #TODO: Clear abbrev out
        if not validate_only:
            with open(restart_json_path,'w') as json_out:
                out_dict = { 'abbrev_to_IUPAC': abbrev_to_IUPAC,
                        'abbrev_to_type': abbrev_to_type,
                        'abbrev_to_MW': dictFix(abbrev_to_MW),
                        'abbrev_to_NopH': abbrev_to_NopH,
                        'abbrev_to_polSurf': abbrev_to_polSurf,
                        'abbrev_to_msacc': abbrev_to_msacc,
                        'abbrev_to_msdon': abbrev_to_msdon,
                        'abbrev_to_nonsaltIUPAC': abbrev_to_nonsaltIUPAC,
                        'abbrev_to_ProjectionArea': abbrev_to_ProjectionArea,
                        'skip': skip}
                json.dump(out_dict,json_out)
            if test_cg:
                sys.exit("completed test of compound guide file")
        else:
            with open(smiles_json_path,'w') as smiles_out:
                json.dump(IUPAC_to_SMILES,smiles_out) 
    ml_convert = json.load(open(mlConvert_json_path))
    with open(rxn_database_path, 'r') as rxn_file:
        rxn_reader = csv.reader(rxn_file)
        rxn_reader.next()
        if not validate_only:
            out_file = open(expanded_csv_path, 'wb')
            out_writer = csv.writer(out_file)
            hdrs = ['XXXtitle', 'XXXinorg1', 'XXXinorg1mass', 
            'XXXinorg1moles', 'XXXinorg2', 'XXXinorg2mass', 
            'XXXinorg2moles', 'XXXinorg3', 'XXXinorg3mass','XXXinorg3moles', 'XXXorg1', 'XXXorg1mass',
            'XXXorg1moles', 'XXXorg2', 'XXXorg2mass', 
            'XXXorg2moles', 'XXXoxlike1', 'XXXoxlike1mass', 
            'XXXoxlike1moles', 'Temp_max', 'time', 'slowCool', 'pH', 
            'leak', 'numberInorg', 'numberOrg', 'numberOxlike', 
            'numberComponents', 'inorgavgpolMax', 
            'inorgrefractivityMax', 'inorgmaximalprojectionareaMax', 
            'inorgmaximalprojectionradiusMax', 
            'inorgmaximalprojectionsizeMax', 
            'inorgminimalprojectionareaMax', 
            'inorgminimalprojectionradiusMax', 
            'inorgminimalprojectionsizeMax', 
            'inorgavgpol_pHdependentMax', 
            'inorgmolpolMax', 'inorgvanderwaalsMax', 'inorgASAMax', 
            'inorgASA+Max', 'inorgASA-Max', 'inorgASA_HMax', 
            'inorgASA_PMax', 'inorgpolarsurfaceareaMax', 
            'inorghbdamsaccMax', 'inorghbdamsdonMax', 
            'inorgavgpolMin', 'inorgrefractivityMin', 
            'inorgmaximalprojectionareaMin', 
            'inorgmaximalprojectionradiusMin', 
            'inorgmaximalprojectionsizeMin', 
            'inorgminimalprojectionareaMin', 
            'inorgminimalprojectionradiusMin', 
            'inorgminimalprojectionsizeMin', 
            'inorgavgpol_pHdependentMin', 
            'inorgmolpolMin', 'inorgvanderwaalsMin', 'inorgASAMin', 
            'inorgASA+Min', 'inorgASA-Min', 'inorgASA_HMin', 
            'inorgASA_PMin', 'inorgpolarsurfaceareaMin', 
            'inorghbdamsaccMin', 'inorghbdamsdonMin', 
            'inorgavgpolArithAvg', 'inorgrefractivityArithAvg', 
            'inorgmaximalprojectionareaArithAvg', 
            'inorgmaximalprojectionradiusArithAvg', 
            'inorgmaximalprojectionsizeArithAvg', 
            'inorgminimalprojectionareaArithAvg', 
            'inorgminimalprojectionradiusArithAvg', 
            'inorgminimalprojectionsizeArithAvg', 
            'inorgavgpol_pHdependentArithAvg', 'inorgmolpolArithAvg', 
            'inorgvanderwaalsArithAvg', 'inorgASAArithAvg', 
            'inorgASA+ArithAvg', 'inorgASA-ArithAvg', 
            'inorgASA_HArithAvg', 'inorgASA_PArithAvg', 
            'inorgpolarsurfaceareaArithAvg', 
            'inorghbdamsaccArithAvg', 'inorghbdamsdonArithAvg', 
            'inorgavgpolGeomAvg', 'inorgrefractivityGeomAvg', 
            'inorgmaximalprojectionareaGeomAvg', 
            'inorgmaximalprojectionradiusGeomAvg', 
            'inorgmaximalprojectionsizeGeomAvg', 
            'inorgminimalprojectionareaGeomAvg', 
            'inorgminimalprojectionradiusGeomAvg', 
            'inorgminimalprojectionsizeGeomAvg', 
            'inorgavgpol_pHdependentGeomAvg', 
            'inorgmolpolGeomAvg', 'inorgvanderwaalsGeomAvg', 
            'inorgASAGeomAvg', 'inorgASA+GeomAvg', 
            'inorgASA-GeomAvg', 'inorgASA_HGeomAvg', 
            'inorgASA_PGeomAvg', 'inorgpolarsurfaceareaGeomAvg', 
            'inorghbdamsaccGeomAvg', 'inorghbdamsdonGeomAvg', 
            'orgavgpolMax', 'orgrefractivityMax', 
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
            'inorgacc-waterdonratio', 'inorgdon-wateraccratio',
            'org-water-moleratio', 'orgacc-waterdonratio',
            'orgdon-wateraccratio', 'inorg-org-moleratio',
            'inorgacc-orgdonratio', 'inorgdon-orgaccratio',
            'notwater-water-moleratio', 'notwateracc-waterdonratio', 
            'notwaterdon-wateraccratio', 'purity', 'outcome']
            out_writer.writerow(hdrs)
            cnt = 0
            for rxn in rxn_reader:
                if cnt % 100000 == 0: ###Necessary or decrease "counter"?
                    print "---- Played with {} reactions so far...".format(cnt)
                cnt += 1
                compound = ["x","x","x","x","x"]
                mass = ["-1","-1","-1","-1","-1"]
                unit = ["","","","",""]
                (title, compound[0], mass[0], unit[0], compound[1], mass[1], unit[1],
                        compound[2], mass[2], unit[2], compound[3], mass[3],
                        unit[3], compound[4], mass[4], unit[4], Tmax,
                        time,  pH, slowCool, leak, outcome, purity) = rxn[:23]
 
                mass = map(lambda x: x[0] + x[1], zip(mass, unit))
                try:
                    compound = [''.join(s.split()) for s in compound]
                    mass = [s.replace('x','-1').replace('z','0') for s in mass]
                    for i in range(5):
                        if compound[i] in skip.keys():
                            compound[i] = 'x'
                            mass[i] = "-1"
                        elif "mL" in mass[i]:
                            assert(compound[i] in ml_convert.keys())
                            mass[i] = float(fixmL(mass[i]))*ml_convert[compound[i]] #TODO
                        elif mass[i] == '' or mass[i] == 'g':
                            mass[i] = "-1"
                        elif mass[i][-1] == 'g':
                            mass[i] = mass[i][:-1]
                        if compound[i] == "":
                            compound[i] = 'x'
                    mass = [float(m) for m in mass] 
                    Tmax = float(Tmax)
                    time = float(time)
                    try:
                        pH = float(pH)
                        if not (pH >= 0 and pH < 14):
                            raise Exception("pH out of bounds")
                    except Exception as e:
                        sys.stderr.write("WARNING: failed pH conversion: {}\n".format(pH))
                        continue
                    if not purity:
                        purity = str(1)
                    purity = int(purity)
                    comb = False 
                    if not comb:
                        if not outcome:
                            outcome = 1
                        outcome = int(outcome)
                    else:
                        if not outcome:
                            outcome = 1
                        outcome = int(str(outcome) + str(purity))
                except Exception as e:
                    exc_tb = sys.exc_info()[2]
                    sys.stderr.write("WARNING: failed to construct rxn. %s %s %s, line number: %s\n" % (str(e),str(type(e)),str(rxn),str(exc_tb.tb_lineno)))
                    continue
                output = [title]
                if not outcome:
                    outcome = 1
                isWater = 0
                keepList = [1,1,1,1,1]
                organicList = [0,0,0,0,0]
                inorganicList = [0,0,0,0,0]
                waterList = [0,0,0,0,0]
                oxalateList = [0,0,0,0,0]
                (nInorg, nOrg, nOxlike) = (0,0,0)
                try:
                    for i in range(5):
                        compound[i] = ''.join(compound[i].split())
                        if compound[i].lower() == "water":
                            isWater = 1
                        elif compound[i].lower() == 'x' or compound[i] in skip.keys():
                            keepList[i] = 0
                        comp_IUPAC = compound[i] in abbrev_to_IUPAC.keys()
                        if keepList[i] and not comp_IUPAC:
                            if not compound[i].replace('-','') in abbrev_to_IUPAC.keys():
                                sys.stderr.write("WARNING: encounter unknown \
                                        abbreviation %s in entry %s. \
                                        SKIPPING.\n" % (compound[i], title))
                                raise Exception('skip')
                except Exception as e:
                    sys.stderr.write('Exception: %s\n' % str(e))
                    continue
                if validate_only:
                    continue

                for i in range(5):
                    if keepList[i]:
                        if abbrev_to_type[compound[i]] == 'i':
                            inorganicList[i] = 1
                            nInorg += 1
                        elif abbrev_to_type[compound[i]] == 'o':
                            organicList[i] = 1
                            nOrg += 1
                        elif abbrev_to_type[compound[i]] == 'ox':
                            oxalateList[i] += 1
                            nOxlike += 1
                        elif abbrev_to_type[compound[i]] == 'w':
                            waterList[i] = 1
                            isWater = i
                if not isWater:
                    sys.stderr.write('WARNING: no water present for \
                            reaction %s\n' % title)
                    continue
                if pH  == 'x':
                    sys.stderr.write("WARNING: encountered unknown pH in entry \
                            %s while parsing %s. SKIPPING THIS ENTRY\n \
                            " % (title, rxn_database_path))
                    continue
                pH = math.ceil(float(pH))
                if pH < 1:
                    sys.stderr.write("Merp!")
                    continue
                bar = ''
                compoundAcc = [-1,-1,-1,-1,-1]
                compoundDon = [-1,-1,-1,-1,-1]
                compoundMoles = [-1,-1,-1,-1,-1]

                orgDetails = []
                inorgDetails = []
                oxlikeDetails = []

                for i in range(5):
                    if not keepList[i]:
                        continue
                    a = abbrev_to_MW[compound[i]]
                    assert(a != 0)
                    compoundMoles[i] = mass[i]/a
                    if isWater and not compoundMoles[i]:
                        sys.stderr.write("Warning: water moles zero, Skipping.\n")
                    r = [compound[i],mass[i],compoundMoles[i]]
                    if inorganicList[i]:
                        inorgDetails += r
                    if organicList[i]:
                        orgDetails += r
                    if oxalateList[i]:
                        oxlikeDetails += r
                for i in [0,1]: #hack, now broken
                    if len(inorgDetails) + 1 < 9:
                        inorgDetails += [-1,-1,-1]
                    if len(orgDetails) + 1 < 6:
                        orgDetails += [-1,-1,-1]
                    if len(oxlikeDetails) + 1 < 3:
                        oxlikeDetails += [-1,-1,-1]

                output += inorgDetails + orgDetails + oxlikeDetails
                output += [Tmax, time, slowCool.lower(), pH, leak.lower(), nInorg,
                        nOrg, nOxlike, nInorg+nOrg+nOxlike]
                compoundProperties = [[] for k in range(5)]
                for j in range(5):
                    if not keepList[j]:
                        compoundProperties[j] = [-1 for k in range(19)]
                        continue
                    a = compound[j]
                    compoundProperties[j] += abbrev_to_NopH[compound[j]]
                    compoundProperties[j] += abbrev_to_ProjectionArea[compound[j]]
                    bar = abbrev_to_polSurf[compound[j]]
                    compoundProperties[j] += [bar[k] for k in range(9*int(pH) -9,
                                9*int(pH))]
                    bar = abbrev_to_msacc[a]
                    compoundAcc[j] = bar[int(pH)-1]
                    if abs(compoundAcc[j]) < 0.01:
                        compoundAcc[j] = 0.0001 # to prevent x/0
                    bar = abbrev_to_msdon[a]
                    compoundDon[j] = bar[int(pH)-1]
                    if abs(compoundDon[j]) < 0.01:
                        compoundDon[i] = 0.0001 #ibid
                    compoundProperties[j] += [compoundAcc[j], compoundDon[j]]
                if nInorg > 0:
                    output += distList(inorganicList, compoundProperties) #questionable. 387 
                else:
                    output += [-1 for i in range(76)] 
                if nOrg > 0:
                    output += distList(organicList, compoundProperties)
                else:
                    output += [-1 for i in range(76)]
                if nOxlike>0:
                    output += distList(oxalateList, compoundProperties)
                else:
                    output += [-1 for i in range(76)]

                #inorganic/water mole ratio
                totalInorganic = 0.0
                for k in range(5):
                    totalInorganic += compoundMoles[k]*inorganicList[k]
                totalWater = compoundMoles[isWater]
                if totalWater == 0:
                    totalWater = 0.00001
                output.append(totalInorganic/totalWater)

                #inorganic/water accepter-on-inorganic ratio
                totalInorganic = 0.0
                for k in range(5):
                    totalInorganic += compoundMoles[k]*inorganicList[k]*compoundAcc[k]
                if not compoundMoles[isWater]:
                    compoundMoles[isWater] = 0.00001 
                if not compoundDon[isWater]:
                    compoundDon[isWater] = 0.00001 
                if not compoundAcc[isWater]:
                    compoundAcc[isWater] = 0.00001 
                output.append(totalInorganic/(compoundMoles[isWater]*compoundDon[isWater]))

                #inorganic/water donor-on-inorganic ratio
                totalInorganic = 0.0
                for k in range(5):
                    totalInorganic += compoundMoles[k]*inorganicList[k]*compoundDon[k]
                output.append(totalInorganic/(compoundMoles[isWater]*compoundAcc[isWater]))

                # organic / water mole ratio
                totalOrganic = 0.0
                for k in range(5):
                    totalOrganic += compoundMoles[k]*organicList[k]
                output.append(totalOrganic/totalWater)

                #organic/water acceptor-on-inorganic ratio
                totalOrganic = 0.0
                for k in range(5):
                    totalOrganic += compoundMoles[k]*organicList[k]*compoundAcc[k]
                output.append(totalOrganic/(compoundMoles[isWater]*compoundDon[isWater]))

                #organic/water donor-on-inorganic ratio
                totalOrganic = 0.0
                for k in range(5):
                    totalOrganic += compoundMoles[k]*organicList[k]*compoundDon[k]
                output.append(totalOrganic/(compoundMoles[isWater]*compoundAcc[isWater]))

                #inorrganic/organic mole-ratio
                totalInorganic = 0.0
                totalOrganic = 0.0
                for k in range(5):
                    totalOrganic += compoundMoles[k]*organicList[k]
                    totalInorganic += compoundMoles[k]*inorganicList[k]
                if not totalOrganic:
                    totalOrganic = 0.00001 
                output.append(totalInorganic/totalOrganic)

                #inorganic/organic acceptor-on-inorganic ratio
                totalInorganic = 0.0
                totalOrganic = 0.0
                for k in range(5):
                    totalOrganic += compoundMoles[k]*organicList[k]*compoundDon[k]
                    totalInorganic += compoundMoles[k]*inorganicList[k]*compoundAcc[k]
                if totalOrganic == 0.0:
                    totalOrganic = 0.00001
                output.append(totalInorganic/totalOrganic)

                #inorganic/organic donor-on-inorganic ratio
                totalInorganic = 0.0
                totalOrganic = 0.0
                for k in range(5):
                    totalOrganic += compoundMoles[k]*organicList[k]*compoundAcc[k]
                    totalInorganic += compoundMoles[k]*inorganicList[k]*compoundDon[k]
                if not totalOrganic:
                    totalOrganic = 0.00001
                output.append(totalInorganic/totalOrganic)
                
                #not water/water mole ratio
                totalInorganic = 0.0
                totalOrganic = 0.0
                for k in range(5):
                    totalOrganic += compoundMoles[k]*organicList[k]
                    totalInorganic += compoundMoles[k]*inorganicList[k]
                if not totalWater:
                    totalWater = 0.00001
                output.append((totalInorganic+totalOrganic)/totalWater)

                #not water/water acceptor on not-water
                totalInorganic = 0.0
                totalOrganic = 0.0
                for k in range(5):
                    totalOrganic += compoundMoles[k]*organicList[k]*compoundAcc[k]
                    totalInorganic += compoundMoles[k]*inorganicList[k]*compoundAcc[k]
                if not compoundMoles[isWater]:
                    compoundMoles[isWater] = 0.00001
                if not compoundDon[isWater]:
                    compoundDon[isWater] = 0.00001
                output.append((totalOrganic+totalInorganic)/(compoundMoles[isWater]*compoundDon[isWater]))

                #not water/water donor-on-not-water ratio
                totalInorganic = 0.0
                totalOrganic = 0.0
                for k in range(5):
                    totalOrganic += compoundMoles[k]*organicList[k]*compoundDon[k]
                    totalInorganic += compoundMoles[k]*inorganicList[k]*compoundDon[k]
                if compoundMoles[isWater] == 0.0:
                    compoundMoles[isWater] = 0.00001
                if compoundAcc[isWater] == 0.0:
                    compoundAcc[isWater] = 0.00001
                output.append((totalOrganic+totalInorganic)/(compoundMoles[isWater]*compoundAcc[isWater]))

                output.append(purity)
                output.append(outcome)
                out_writer.writerow(map(lambda s: "{0:.4f}".format(s) if type(s) == float else s,  output))
    if validate_only:
        print "----    Validation complete!"
    else:
        print "---- constructDescriptorTable.py: Success!"
    


