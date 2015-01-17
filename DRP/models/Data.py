from django.db import models
from DRP.data_config import CONFIG

from Lab_Group import Lab_Group
from DataCalc import DataCalc
from django.contrib.auth.models import User

#Many data are saved per lab group. Each data represents one submission.
class Data(models.Model):
  class Meta:
    app_label = "DRP"

  ref = models.CharField("Reference", max_length=12)

  #List Fields
  for i in CONFIG.reactant_range():
    exec("reactant_{0} = models.CharField(\"Reactant {0}\", max_length=30)".format(i))
    exec("quantity_{0} = models.CharField(\"Quantity {0}\", max_length=10)".format(i))
    exec("unit_{0} = models.CharField(\"Unit {0}\", max_length=4)".format(i))

  temp = models.CharField("Temperature", max_length=10)
  time = models.CharField("Time", max_length=10) ###
  pH = models.CharField("pH", max_length=5)

  #Yes/No/? Fields:
  slow_cool = models.CharField("Slow Cool", max_length=10)
  leak = models.CharField("Leak", max_length=10)
  outcome = models.CharField("Outcome", max_length=1)
  purity = models.CharField("Purity", max_length=1)

  notes = models.CharField("Notes", max_length=200, blank=True)

  #Self-assigning Fields:
  calculations = models.ForeignKey(DataCalc, unique=False, blank=True, null=True)
  calculated_pH = models.BooleanField(default=False)
  calculated_temp = models.BooleanField(default=False)
  calculated_time = models.BooleanField(default=False)

  atoms = models.CharField("Atoms", max_length=30, blank=True)

  user = models.ForeignKey(User, unique=False)
  lab_group = models.ForeignKey(Lab_Group, unique=False)
  creation_time_dt = models.DateTimeField("Created", null=True, blank=True)
  is_valid = models.BooleanField("Valid", default=False)

  #Categorizing Fields:
  public = models.BooleanField("Public", default=False)
  duplicate_of = models.CharField("Duplicate", max_length=12, null=True, blank=True)
  recommended = models.CharField("Recommended", max_length=10)

  def __unicode__(self):
    return u"{} -- (LAB: {})".format(self.ref, self.lab_group.lab_title)


  def get_calculations_dict(self, include_lab_info=False, force_recalculate=False,
                            preloaded_cg=None, debug=False,
                            preloaded_abbrev_map=None):
    from model_building.load_data import create_expanded_datum_field_list
    from model_building.rxn_calculator import headers

    try:
      if not self.calculations or force_recalculate:
        # Create the extended calculations.
        calcList = create_expanded_datum_field_list(self, preloaded_cg=preloaded_cg,
                                                    preloaded_abbrev_map=preloaded_abbrev_map)

        calcDict = {key:make_float(val) for key,val in zip(headers, calcList)}

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

    except Exception as e:
      # If a field isn't present in the calcDict, update the calculation.
      calcDict = self.get_calculations_dict(include_lab_info=include_lab_info,
                                            force_recalculate=True,
                                            preloaded_cg=preloaded_cg,
                                            preloaded_abbrev_map=preloaded_abbrev_map)
      return [calcDict[field] for field in headers]


  def to_list(self):
    all_fields = get_model_field_names(model="Data", collect_ignored = True)
    fields_to_exclude = {"lab_group", "atoms"}
    headings = [field for field in all_fields if field not in fields_to_exclude]

    return [getattr(datum,field) for field in headings]



def get_lab_Data(lab_group):
  from DRP.models import get_Lab_Group
  lab_group = get_Lab_Group(lab_group)
  return Data.objects.filter(lab_group=lab_group).order_by("creation_time_dt")


def get_ref_set(lab_group, reset_cache=True):
  from DRP.cacheFunctions import get_cache, set_cache

  ref_set = get_cache(lab_group, "DATAREFS")
  if not ref_set or reset_cache:
    ref_set = set(get_lab_Data(lab_group).values_list('ref', flat=True))
    set_cache(lab_group, "DATAREFS", ref_set)

  return ref_set
