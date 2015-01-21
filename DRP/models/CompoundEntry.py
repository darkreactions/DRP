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


  def get_atoms(self):
    import rdkit.Chem as Chem

    if self.custom and not self.smiles:
      message = "Cannot get SMILES of custom compound: '{}'".format(self.abbrev)
      raise Exception(message)
    elif not self.smiles:
      message = "Compound has no SMILES: '{}'".format(self.abbrev)
      raise Exception(message)


    mols = Chem.MolFromSmiles(self.smiles,sanitize=False)
    if mols == None:
      return []

    #TODO: Incorporate hydrogens into model.
    #if show_hydrogen:
    # try:
    #  mols = Chem.AddHs(mols) ###PRECONDITION?
    # except:
    #  pass

    return [atom.GetSymbol() for atom in mols.GetAtoms()]



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


