
import sys, os
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.settings import BASE_DIR
import load_cg,json


def load(lab_group=None, with_headings=True):
        from DRP.models import get_good_rxns
        return get_good_rxns(lab_group=lab_group, with_headings=with_headings)

#Translate the abbrevs to the full compound names.
def fix_abbrevs(rxnList, abbrev_map=None):
  # Create the "abbrev_map" which will translate
  if not abbrev_map:
    abbrev_map = get_abbrev_map()

  indexes = [1,4,7,10,13]
  for i in indexes:
    if rxnList[i] in abbrev_map:
      rxnList[i] = abbrev_map[rxnList[i]]
  return rxnList


def get_abbrev_map():
  from DRP.models import CompoundEntry as c

  entries = c.objects.all()
  abbrev_map = dict()
  compound_set = set()

  for e in entries:
    abbrev_map[e.abbrev] = e.compound
    compound_set.add(e.compound)

  abbrev_map[''] = ''
  return abbrev_map


def get_feature_vectors(lab_group=None, cg = None, ml_convert = None, keys = None):
        raw = load(lab_group)
        return convert_to_feature_vectors(raw,cg, ml_convert, keys = keys)



def create_expanded_datum_field_list(datum, preloaded_cg=None,
                                            preloaded_abbrev_map=None):
  from DRP.model_building.parse_rxn import parse_rxn

  #Convert the datum into a datumList (ie: a list of field values).
  dirtyDatumList = datum.to_list(keep_foreign_keys=False)

  #Ignore the "ref" field
  datumList = fix_abbrevs( dirtyDatumList, abbrev_map=preloaded_abbrev_map )

  #Grab the Compound Guide and mL Conversion Guide for parse_rxn.
  compoundGuide = load_cg.get_cg() if not preloaded_cg else preloaded_cg
  ml_convert = json.load(open("{}/DRP/model_building/mlConvert.json".format(BASE_DIR)))

  #Actually parse the reaction and spit out an "expanded" one.
  calculations = parse_rxn(datumList, compoundGuide, ml_convert)

  return calculations


def convert_to_feature_vectors(raw, cg = None, ml_convert = None, keys = None):
  if not cg:
    cg = load_cg.get_cg()
  if not ml_convert:
    ml_convert = json.load(open("{}/DRP/model_building/mlConvert.json".format(BASE_DIR)))
  import parse_rxn

  transformed = []
  failed = 0
  keys = []
  for row in raw:
    try:
      calculations = parse_rxn.parse_rxn(row, cg, ml_convert)
      transformed.append(calculations)
      keys.append(create_key(row))
    except Exception as e:
      failed += 1
      print "ERROR convert_to_feature_vectors: {}".format(e)
  if failed:
    print "{0} failed out of {1} total".format(failed, len(raw))
  remove_XXX(transformed)

  #for r in transformed:
  #  del r[-2]

  if keys:
    return transformed, keys

  return transformed


def create_reactant_keys(data, headers):
  from DRP.model_building.rxn_calculator import reactant_fields

  # Variable Setup
  blacklist = {"water", "", -1}
  keys = []

  reactant_indexes = [i for i, h in enumerate(headers) if h in reactant_fields]
  if not reactant_indexes:
    raise Exception("No `reactant_indexes` specified. Splitting aborted.")

  for row in data:
    reactants = [row[i] for i in reactant_indexes if row[i] not in blacklist]
    keys.append( tuple(sorted(reactants)) )

  return keys


def create_key(line):
        #print line
        key = [line[0], line[3], line[6], line[9], line[12]]
        key = [r for r in key if r.lower() != 'water' and r != ""]
        key.sort()
        return tuple(key)


def get_feature_vectors_by_triple(lab_group=None, cg = None, ml_convert = None):
        import parse_rxn
        from DRP.models import CompoundEntry
        if not cg:
                cg = load_cg.get_cg()
        if not ml_convert:
                ml_convert = json.load(open("{}/DRP/model_building/mlConvert.json".format(BASE_DIR)))

        raw = load(lab_group)[1:] #daniel (taking out the headers)
        def rect(arr_2d): #daniel
            return all(len(i) == len(arr_2d[0]) for i in arr_2d) #daniel
        import sys, os #daniel
        sys.stdout.write("all rectangular?: " + str(all(rect(item) for item in raw)) + "\n") #daniel
        os.system('echo "' + "are they rectangular?" + '"|espeak') #daniel
        raise(Exception("are they rectangular?")) #daniel
        transformed = []

        triple_to_rxn_list = dict()
        import sys, os #daniel
        blahh = "here in get_feature_vectors_by_triple in DRP model_building load_data dot py" #daniel
        sys.stdout.write(blahh + "\n") #daniel
        os.system('echo "' + blahh + '"|espeak') #daniel
        sys.stdout.flush() #daniel
        for rxn in raw:
                try:
                        triple = rxn_to_triple(rxn, cg)
                except Exception as e:
                        continue
                if triple not in triple_to_rxn_list:
                        triple_to_rxn_list[triple] = []
                triple_to_rxn_list[triple].append(rxn)
        for triple in triple_to_rxn_list:
                rxn_list = triple_to_rxn_list[triple]
                transformed = []
                for row in rxn_list:
                        try:
                                row_w_cmpd_strs = [x.compound if type(x) == CompoundEntry else x for x in row] #daniel
                                transformed.append(parse_rxn.parse_rxn(row_w_cmpd_strs, cg, ml_convert))
                        except Exception as e:
                                print "EXCEPTION: {}".format(e)

                remove_XXX(transformed)
                triple_to_rxn_list[triple] = transformed

        return triple_to_rxn_list

def collapse_triples(dataset):
        unknown = []
        del_triples = []
        for triple in dataset:
                if len(dataset[triple]) < 6:
                        unknown += dataset[triple]
                        del_triples.append(triple)
        for triple in del_triples:
                del dataset[triple]

        dataset['unknown'] = unknown


# Daniel: this seems originally written to take abbrevs (CompoundEntry.abbrev, abbreviation string for a compound), 
#   but is used only one place, and is now taking CompoundEntry objects. Rewriting to expect CompoundEntry objects instead.
def rxn_to_triple(rxn, cg):
        from DRP.compoundGuideFunctions import translate_reactants
        r = rxn
        #import sys #daniel
        #sys.stdout.write("rxn: " + str(r) + "\n") #daniel
        compounds = filter(lambda x: x != None and x.compound.lower() != 'water' and x.compound.lower() != '', [r[1], r[4], r[7], r[10], r[13]]) #daniel
        #sys.stdout.write("compounds: " + str(compounds) + "\n") # daniel
        #sys.stdout.write("compounds type: " + str(type(compounds[0])) + "\n") #daniel

        #TODO: Should pass the specific LabGroup into here!
        # Daniel: this appears to be pointless now that it takes CompoundEntries?
        #compounds = translate_reactants("Norquist Lab", compounds, #daniel
        #                                onlyAbbrevs=True, #daniel
        #                                direction="abbrev to compound") #daniel
        #sys.stdout.write("compounds after translate: " + str(compounds) + "\n") #daniel
        #sys.stdout.write("compounds type after translate: " + str(type(compounds[0])) + "\n") #daniel
        #sys.stdout.flush()
        #import sys, os #daniel
        #sys.stdout.write("translate_reactants completes successfully" + "\n") #daniel
        for compound in compounds:
                if compound.compound.lower() not in cg: #daniel; added .compound.lower()
                        raise Exception("Unknown compound: {0}".format(compound))
        return tuple(sorted(compounds))



def remove_XXX(row):
  import rxn_calculator
  dist = 0
  end = rxn_calculator.headers.index('outcome') + 1
  for hdr in rxn_calculator.headers:
    if "XXX" in hdr:
      dist += 1
  row = row[dist:end]
  return row
