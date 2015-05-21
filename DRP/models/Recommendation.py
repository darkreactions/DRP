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

  #Giving Rec objects a ref field
  ref = models.CharField("Reference", max_length=30, unique=False, default="")

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

  def to_list(self, keep_foreign_keys=True):
    from methods import get_model_field_names

    # Variable Setup
    fields_to_exclude = {"lab_group", "atoms"}

    # Get the relevant headings.
    all_fields = get_model_field_names(model="Data",
                                       collect_ignored = True)
    headings = [field for field in all_fields if field not in fields_to_exclude]

    row = [getattr(self,field) for field in headings if hasattr(self,field)]

    # Convert CompoundEntry ForeignKeys to their respective "compounds".
    if not keep_foreign_keys:
      row = [elem if type(elem)!=CompoundEntry else elem.compound for elem in row]
      row = map(lambda elem: elem if elem is not None else "", row)

    return row

  def get_calculations_dict(self, include_lab_info=False, force_recalculate=False,
                            preloaded_cg=None, debug=False,
                            preloaded_abbrev_map=None):

    def _make_float(raw):
      # Convert any number-like strings to floats.
      try:
        return float(raw)
      except:
        return raw


    from DRP.model_building.load_data import create_expanded_datum_field_list
    from DRP.model_building.rxn_calculator import headers
    import json

    try:
      if not self.calculations or force_recalculate:
        # Create the extended calculations.
        calcList = create_expanded_datum_field_list(self, preloaded_cg=preloaded_cg,
                                                    preloaded_abbrev_map=preloaded_abbrev_map)

        calcDict = {key:_make_float(val) for key,val in zip(headers, calcList)}

        # Prepare a new DataCalc object.
        newDataCalc = DataCalc(contents=json.dumps(calcDict))
        newDataCalc.save()

        # Create the ForeignKey between the new DataCalc and this Datum.
        self.calculations = newDataCalc
        self.is_valid = True
        self.save()

        final_dict = calcDict

      else:
        # Load the result from the database if it is already present.
        final_dict = self.calculations.make_json()


      if include_lab_info:
        final_dict.update({
                          "lab_title":self.lab_group.lab_title,
                          "creation_time_dt":str(self.creation_time_dt),
                          })

      # Make sure all of the keys are present.
      missing_keys = not set(headers).issubset(set(final_dict.keys()))

      if type(final_dict)!=dict or missing_keys:
        # If the final_dict is in the wrong format, recalculate it.
        return self.get_calculations_dict(include_lab_info=include_lab_info,
                                          force_recalculate=True)

      if self.nonsense and "outcome" in final_dict:
        final_dict["outcome"] = 0

      return final_dict

    except Exception as e:
      self.is_valid = False
      self.save()
      raise Exception("(get_calculations_dict) {}".format(e))

  def get_calculations_list(self, include_lab_info=False, preloaded_cg=None,
                           preloaded_abbrev_map=None):
    from DRP.model_building.rxn_calculator import headers

    if include_lab_info:
      headers.extend(["lab_title", "creation_time_dt"])

    try:
      calcDict = self.get_calculations_dict(include_lab_info=include_lab_info,
                                            preloaded_cg=preloaded_cg,
                                            preloaded_abbrev_map=preloaded_abbrev_map)
      return [calcDict[field] for field in headers]

    except:
      # If a field isn't present in the calcDict, update the calculation.
      calcDict = self.get_calculations_dict(include_lab_info=include_lab_info,
                                            force_recalculate=True,
                                            preloaded_cg=preloaded_cg,
                                            preloaded_abbrev_map=preloaded_abbrev_map)
      return [calcDict[field] for field in headers]


def get_recommendations(lab_group):
  return Recommendation.objects.filter(lab_group=lab_group)


def get_recommendations_by_date(lab_group):
  return get_recommendations(lab_group).order_by("-date_dt")

