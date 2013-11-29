####################
################### This File is not Stable.
####################
import rdkit.Chem as Chem
import chemspipy
from data_config import CONFIG

def chemspider_lookup(cg_entry):
 chemspi_query = ""
 try:
  #Accept either a CompoundEntry object or a dict with the valid fields
  search_fields = ["CAS_ID", "compound"]
  if type(cg_entry)==dict:
   query_criteria = [cg_entry.get(i) for i in search_fields if cg_entry.get(i)]
  elif type(cg_entry)==CompoundEntry:
   query_criteria = [getattr(cg_entry, i) for i in search_fields]
  elif type(cg_entry)==str:
   query_criteria = [cg_entry]

  #Works for both dicts and CG_entries
  assert(query_criteria)
 except:
  raise Exception("chemspider_lookup accepts strings, dicts, or CompoundEntries")

 for i in query_criteria:
  if i: #ChemSpider doesn't like empty requests.
   chemspi_query = chemspipy.find_one(i)
  if chemspi_query: break

 if chemspi_query:
  return chemspi_query.imageurl, chemspi_query.smiles, chemspi_query.molecularweight
 else:
  raise Exception("Could not find compound on ChemSpider!")

def get_atoms_from_compound(CG_entry = None):
 return get_atoms_from_smiles(CG_entry.smiles)

def get_atoms_from_smiles(smiles, show_hydrogen=False):
 if not smiles:
  print "No smiles found!"
  return []

 mols = Chem.MolFromSmiles(str(smiles),sanitize=False)
 if mols == None:
  return []
 #TODO: Incorporate hydrogens into model.
 #if show_hydrogen:
 # try:
 #  mols = Chem.AddHs(mols) ###PRECONDITION?
 # except:
 #  pass
 atoms = mols.GetAtoms()
 return [atom.GetSymbol() for atom in atoms]

def collect_CGs_by_abbrevs(lab_group, abbrev_list):
 CG_list = []
 for i in abbrev_list:
  query = CompoundEntry.objects.filter(lab_group=lab_group, abbrev=i)
  if query.exists(): #Don't append index an empty query; Django gets annoyed.
   CG_list.append(query[0])
 return CG_list

def get_smiles_from_CG_list(CG_list):
 smiles_list = [i.smiles for i in CG_list]
 return smiles_list

def condense_smiles_list_to_atoms(smiles_list, show_hydrogen = False):
 atoms_list = []
 for i in smiles_list:
  atoms_list += get_atoms_from_smiles(i, show_hydrogen)
 return set(atoms_list)

def get_abbrevs_from_reaction(reaction):
 abbrevs_list = [getattr(reaction, "reactant_{}".format(i)) for i in CONFIG.reactant_range() if getattr(reaction, "reactant_{}".format(i))]
 return abbrevs_list

def get_atom_set_from_abbrevs(lab_group, abbrev_list):
 return condense_smiles_list_to_atoms(
  get_smiles_from_CG_list(
   collect_CGs_by_abbrevs(lab_group, abbrev_list)
   ), show_hydrogen=True
  )

def get_atom_set_from_reaction(reaction):
 return get_atom_set_from_abbrevs(reaction.lab_group, get_abbrevs_from_reaction(reaction))

def get_reactions_with_compound(compound):
 Q_string = ""
 abbrev = compound.abbrev
 for i in CONFIG.reactant_range():
  Q_string += "Q(reactant_{}=\"{}\")|".format(i, abbrev)

 if Q_string:
  Q_string = Q_string[:-1] #Remove the last trailing "|" if applicable.

 return eval("Data.objects.filter(lab_group=compound.lab_group).filter({})".format(Q_string))

def update_compound(compound, update_data=True):
 try:
  #Update the CG entry itself. Make sure "inorg" types don't query ChemSpider.
  try:
   compound.image_url, compound.smiles, compound.mw = chemspider_lookup(compound)
   print compound.compound_type, compound.abbrev
  except:
   if compound.compound_type=="Inorg":
    compound.image_url, compound.smiles, compound.mw = "","",""
   else:
    raise Exception("Could not find via ChemSpider!")
  compound.save()

  #Update the individual "atom" records on each reaction.
  if update_data:
   update_reactions(compound)
 except Exception as e:
  print "Could not update {}\n\t{}".format(compound, e)

def update_reactions(compound):
 #Update the individual "atom" records on each reaction.
 changed_reactions = get_reactions_with_compound(compound)
 for reaction in changed_reactions:
  #Store the atoms as a string -- not a set.
  reaction.atoms = "".join(get_atom_set_from_reaction(reaction))
  reaction.save()
