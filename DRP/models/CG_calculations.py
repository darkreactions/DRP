from django.db import models

class CG_calculations(models.Model):
  class Meta:
    app_label = "DRP"

  json_data = models.TextField()
  compound = models.CharField(max_length=200)
  smiles = models.CharField(max_length=255, unique=True)
  json = models.TextField(null=True, default="{}")

  def __unicode__(self):
    return u"{} ({})".format(self.compound, self.smiles)


def create_CG_calcs_if_needed(compound, smiles, compound_type, ):
  from CGCalculator import CGCalculator
  from DRP.data_config import CONFIG
  from DRP.models import CompoundEntry
  import time, json

  #Variable Setup
  jchem_path =  CONFIG.jchem_path
  sdf_path = "tmp"

  compound, smiles, compound_type = str(compound), str(smiles), str(compound_type)

  #Only Organics that have smiles may have calculations.
  if compound_type not in {"Org", "Inorg"} or not smiles:
    return

  #Either return an old CG_calculation or a new one.
  print compound_type

  try:
    cgc = CG_calculations.objects.filter(smiles=smiles)[0]
  except:
    #Calculate properties for the CGEntry
    sdf_filename = str(int(time.time())) + filter(str.isalnum, compound)
    #TODO: Speed this up? This is dreadfully slow.
    props = CGCalculator(compound, sdf_filename, smiles, compound_type, jchem_path, sdf_path).get_properties()
    props = json.dumps(props)
    #Store the actual CG_calculation in the database.
    cgc = CG_calculations(json_data=props, compound=compound, smiles=smiles)
    cgc.save()

    #Set the calculations field in each CompoundEntry.
    CompoundEntry.objects.filter(smiles=smiles).update(calculations=cgc)

  return cgc


