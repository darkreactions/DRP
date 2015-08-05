'''Module containing only the Compound Class'''
from django.db import models
from MolDescriptor import MolDescriptor
from ChemicalClass import ChemicalClass
from LabGroup import LabGroup, get_Lab_Group

class Compound(models.Model):
  '''A class for containing data about Compounds used in chemical reactions.
  The assumption is made that all chemicals used are single-species.
  '''
  
  class Meta:
    app_label = "DRP"

  abbrev = models.CharField("Abbreviation", max_length=100)
  '''A local, often nonstandard abbreviation for a compound'''
  name = models.CharField('Name:', max_length=300)
  '''Normally the IUPAC name of the compound, however this may not be the most parsable name (which is preferable)'''
  CAS_ID = models.CharField("CAS ID", max_length=13, blank=True, default="")
  CHEBI_ID = models.PositiveIntegerField('CHEBI ID')
  '''The CHEBI_ID for the compound. This has been added as a zero-cost attribute which may be used for automating chemical semantics later'''
  ChemicalClass = models.ManyToManyField(ChemicalClass)
  '''The class of the compound- examples include Inorganic Salt'''
  CSID = models.PositiveIntegerField('Chemspider ID')
  '''The chemspider ID for the compound- preferable to the CAS_ID since it is not subject to licensing restrictions'''
  custom = models.BooleanField("Custom", default=False)
  '''This flag denotes whether a compound has been added irrespective of other validation.
  This should be restricted to superusers'''
  INCHI = models.TextField('InCHI key', blank=True, default='')
  '''The Inchi key for a compound- a canonical representation of a molecule which is also unique.'''

  smiles= models.TextField('Smiles', blank=True, default='')
  '''A non-canonical string representation of a molecule which cannot be directly used to test for identity
  but is nevertheless useful for calculating descriptors
  '''

  descriptors = models.ManyToManyField(MolDescriptor, through='MolDescriptorValue')
  '''A link to descriptors which have been calculated for this compound. Values for the descriptors are found
  on the MolDescriptorValue model.
  '''

  def __unicode__(self):
    if self.compound == self.abbrev:
      return u"{} (--> same) (LAB: {})".format(self.abbrev, self.lab_group.lab_title)
    return u"{} --> {} (LAB: {})".format(self.abbrev, self.compound, self.lab_group.lab_title)


  def get_atoms(self, fail_soft=False):
    import rdkit.Chem as Chem

    error = ""

    if self.custom and not self.smiles:
      error = "Cannot get SMILES of custom compound: '{}'".format(self.abbrev)
    elif not self.smiles:
      error = "Compound has no SMILES: '{}'".format(self.abbrev)

    if error and not fail_soft:
      raise Exception(error)

    print "I made it" 
    mols = Chem.MolFromSmiles(str(self.smiles),sanitize=False)
    if mols == None:
      return []

    #TODO: Incorporate hydrogens into model.
    #if show_hydrogen:
    # try:
    #  mols = Chem.AddHs(mols) ###PRECONDITION?
    # except:
    #  pass

    return [atom.GetSymbol() for atom in mols.GetAtoms()]

  def create_CG_calcs_if_needed(self):
    from DRP.CGCalculator import CGCalculator
    from DRP.data_config import CONFIG
    import time, json
    from DRP.fileFunctions import createDirIfNecessary
    #Variable Setup
    jchem_path =  CONFIG.jchem_path
    sdf_path = "tmp"
    createDirIfNecessary(sdf_path) 
    #Only Organics that have smiles may have calculations.
    if self.compound_type not in {"Org", "Inorg"} or not self.smiles:
      return

    #Either return an old CG_calculation or a new one.
    try:
      cgc = CG_calculations.objects.filter(smiles=self.smiles)[0]
    except:
      #Calculate properties for the CGEntry
      sdf_filename = str(int(time.time())) + "".join(filter(str.isalnum, str(self.compound)))
      #TODO: Speed this up? This is dreadfully slow.
      props = CGCalculator(self.compound, sdf_filename, self.smiles, self.compound_type, jchem_path, sdf_path).get_properties()
      props = json.dumps(props)
      #Store the actual CG_calculation in the database.
      cgc = CG_calculations(json_data=props, compound=self.compound, smiles=self.smiles)
      cgc.save()

      #Set the calculations field in each CompoundEntry.
      self.calculations=cgc
      self.save()

    return cgc








def get_compound(abbrev, lab_group):
  """
  Returns the CompoundEntry object with a given `abbrev`.
  """

  compounds = CompoundEntry.objects.all()

  if lab_group:
    lab_group = get_Lab_Group(lab_group)
    compounds = compounds.filter(lab_group=lab_group)

  return compounds.filter(abbrev=abbrev).first()


def compound_exists(abbrev, lab_group=""):
  """
  Returns  `True` if a compound is found and `False` if it is not.
  """

  try:
    comp = get_compound(abbrev, lab_group=lab_group)
    return comp is not None
  except Exception as e:
    print e
    return False


def get_lab_CG(lab_query):
  lab_group = get_Lab_Group(lab_query)
  return CompoundEntry.objects.filter(lab_group=lab_group).order_by("compound")


def collect_CG_name_pairs(lab_group, reset_cache=False):
  from DRP.cacheFunctions import get_cache, set_cache

  pairs = get_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS")
  if not pairs or reset_cache:
    compound_guide = get_lab_CG(lab_group)
    pairs = {entry.abbrev: entry.compound for entry in compound_guide}
    set_cache(lab_group, "COMPOUNDGUIDE|NAMEPAIRS", pairs)

  return pairs


