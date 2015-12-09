from django.db import models

from Lab_Group import Lab_Group, get_Lab_Group
from CG_calculations import CG_calculations

class CompoundEntry(models.Model):
  class Meta:
    app_label = "DRP"

  abbrev = models.CharField("Abbreviation", max_length=100)
  compound = models.CharField("Compound", max_length=100)
  CAS_ID = models.CharField("CAS ID", max_length=13, blank=True, default="")
  compound_type = models.CharField("Type", max_length=10)
  image_url = models.CharField("Image URL", max_length=100, blank=True, default="")
  smiles = models.CharField("SMILES", max_length=255, blank=True, default="")
  mw = models.CharField("Molecular Weight", max_length=20, default="")
  custom = models.BooleanField("Custom", default=False)

  lab_group = models.ForeignKey(Lab_Group, unique=False)
  calculations = models.ForeignKey(CG_calculations, unique=False, null=True, default=None)
  calculations_failed = models.BooleanField(default=False)


  def __unicode__(self):
    if self.compound == self.abbrev:
      return u"{} (--> same) (LAB: {})".format(self.abbrev, self.lab_group.lab_title)
    return u"{} --> {} (LAB: {})".format(self.abbrev, self.compound, self.lab_group.lab_title)


  def get_atoms(self, fail_soft=False):
    import rdkit.Chem as Chem

    error = ""

    if self.compound == '' and self.abbrev == '': # This clause added by Daniel, 08-Jul-15
      return []
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


