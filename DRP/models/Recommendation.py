from django.db import models
from django.contrib.auth.models import User

from DRP.data_config import CONFIG

from Lab_Group import Lab_Group
from ModelStats import ModelStats

from Data import Data
from CompoundEntry import CompoundEntry
from DataCalc import DataCalc

class Recommendation(models.Model):
  class Meta:
    app_label = "DRP"

  #Reactant Fields
  for i in CONFIG.reactant_range():
    exec("reactant_fk_{0} = models.ForeignKey(CompoundEntry, max_length=30, default=None, null=True, blank=True, related_name='rec_reactant_key_{0}')".format(i))

    exec("quantity_{0} = models.CharField(\"Quantity {0}\", max_length=10)".format(i))
    exec("unit_{0} = models.CharField(\"Unit {0}\", max_length=4)".format(i))

  score = models.FloatField("Score")
  temp = models.CharField("Temperature", max_length=10)
  time = models.CharField("Time", max_length=10) ###
  pH = models.CharField("pH", max_length=5)

  #Yes/No/? Fields:
  slow_cool = models.CharField("Slow Cool", max_length=10)
  leak = models.CharField("Leak", max_length=10)
  outcome = models.CharField("Outcome", max_length=1)
  purity = models.CharField("Purity", max_length=1)

  #Self-assigning Fields:
  atoms = models.CharField("Atoms", max_length=30, blank=True)
  lab_group = models.ForeignKey(Lab_Group, unique=False)
  model_version = models.ForeignKey(ModelStats, unique=False)
  user = models.ForeignKey(User, unique=False, null=True, blank=True, default=None, related_name="last_user")
  assigned_user = models.ForeignKey(User, unique=False, null=True, blank=True, default=None, related_name="assigned_user")
  seed = models.ForeignKey(Data, unique=False, null=True, blank=True, default=None)
  date_dt = models.DateTimeField("Created", null=True, blank=True)
  complete = models.BooleanField("Complete", default=False)
  seeded = models.BooleanField("From Seed", default=False)

  #Fields for user feedback.
  saved = models.BooleanField("Saved", default=False)
  nonsense = models.BooleanField("Nonsense", default=False)
  hidden = models.BooleanField("Hidden", default=False)
  notes = models.CharField("Notes", max_length=200, blank=True)


  #Foreign Key field
  calculations = models.ForeignKey(DataCalc, unique=False, blank=True, null=True,

                                   on_delete=models.SET_NULL)
  #Give Recommendation a 'to_list' attribute
  def to_list(self):
    from DRP.retrievalFunctions import get_model_field_names
    all_fields = get_model_field_names(model="Data", collect_ignored = True)
    fields_to_exclude = {"lab_group", "atoms"}
    headings = [field for field in all_fields if field not in fields_to_exclude]

    return [getattr(self,field) for field in headings]

def gather_all_nonsense_recs():
  from DRP.model_building.load_data import create_expanded_datum_field_list
  from DRP.model_building.rxn_calculator import headers
  from DRP.model_building.load_cg import get_cg
  from DRP.model_building.load_data import get_abbrev_map

  nonsense = Recommendation.objects.filter(nonsense=True)
  cg = get_cg()
  abbrev_map = get_abbrev_map()

  for i, elem in enumerate(nonsense):
    nonsenseObjectList = create_expanded_datum_field_list(elem, preloaded_cg=cg,preloaded_abbrev_map=abbrev_map) 
    nonsenseDict = {key:_make_float(val) for key, al in zip(headers, nonsenseObjectListi)} 
    nonsenseDict["outcome"] = 0
    newDataCalc = DataCalc(contents=json.dumps(nonsenseDict))
    newDataCalc.save() 
    elem.calculations = newDataCalc
  return nonsense 

      



def get_recommendations(lab_group):
  return Recommendation.objects.filter(lab_group=lab_group)


def get_recommendations_by_date(lab_group):
  return get_recommendations(lab_group).order_by("-date_dt")

