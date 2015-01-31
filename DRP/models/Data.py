from django.db import models
from DRP.data_config import CONFIG

from Lab_Group import Lab_Group
from DataCalc import DataCalc
from CompoundEntry import CompoundEntry

from django.contrib.auth.models import User

import json

#Many data are saved per lab group. Each data represents one submission.
class Data(models.Model):
  class Meta:
    app_label = "DRP"

  ref = models.CharField("Reference", max_length=30, unique=True)

  ###################################
  #         Reactant Fields         #
  ###################################
  reactant_fk_1 = models.ForeignKey(CompoundEntry, max_length=30, default=None, null=True, blank=True, related_name='reactant_key_1')
  quantity_1 = models.CharField("Quantity 1", max_length=10)
  unit_1 = models.CharField("Unit 1", max_length=4)

  reactant_fk_2 = models.ForeignKey(CompoundEntry, max_length=30, default=None, null=True, blank=True, related_name='reactant_key_2')
  quantity_2 = models.CharField("Quantity 2", max_length=10)
  unit_2 = models.CharField("Unit 2", max_length=4)

  reactant_fk_3 = models.ForeignKey(CompoundEntry, max_length=30, default=None, null=True, blank=True, related_name='reactant_key_3')
  quantity_3 = models.CharField("Quantity 3", max_length=10, blank=True)
  unit_3 = models.CharField("Unit 3", max_length=4, blank=True)

  reactant_fk_4 = models.ForeignKey(CompoundEntry, max_length=30, default=None, null=True, blank=True, related_name='reactant_key_4')
  quantity_4 = models.CharField("Quantity 4", max_length=10, blank=True)
  unit_4 = models.CharField("Unit 4", max_length=4, blank=True)

  reactant_fk_5 = models.ForeignKey(CompoundEntry, max_length=30, default=None, null=True, blank=True, related_name='reactant_key_5')
  quantity_5 = models.CharField("Quantity 5", max_length=10, blank=True)
  unit_5 = models.CharField("Unit 5", max_length=4, blank=True)


  temp = models.CharField("Temperature", max_length=10)
  time = models.CharField("Time", max_length=10)
  pH = models.CharField("pH", max_length=5)

  #Yes/No/? Fields:
  slow_cool = models.CharField("Slow Cool", max_length=10)
  leak = models.CharField("Leak", max_length=10)
  outcome = models.CharField("Outcome", max_length=1)
  purity = models.CharField("Purity", max_length=1)

  notes = models.TextField("Notes", blank=True)

  #Self-assigning Fields:
  calculations = models.ForeignKey(DataCalc, unique=False, blank=True, null=True,
                                   on_delete=models.SET_NULL)
  calculated_pH = models.BooleanField(default=False)
  calculated_temp = models.BooleanField(default=False)
  calculated_time = models.BooleanField(default=False)

  atoms = models.CharField("Atoms", max_length=100, blank=True)

  user = models.ForeignKey(User, unique=False)
  lab_group = models.ForeignKey(Lab_Group, unique=False)
  creation_time_dt = models.DateTimeField("Created", null=True, blank=True)
  is_valid = models.BooleanField("Valid", default=False)

  #Categorizing Fields:
  public = models.BooleanField("Public", default=False)
  duplicate_of = models.CharField("Duplicate", max_length=12, null=True, blank=True)
  recommended = models.CharField("Recommended", max_length=10)

  persistent_homologies = models.TextField("Persistent Homologies", blank=True, default="[]")

  def __unicode__(self):
    return u"{} -- (LAB: {})".format(self.ref, self.lab_group.lab_title)


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


  def to_list(self):
    from DRP.retrievalFunctions import get_model_field_names
    all_fields = get_model_field_names(model="Data", collect_ignored = True)
    fields_to_exclude = {"lab_group", "atoms"}
    headings = [field for field in all_fields if field not in fields_to_exclude]

    return [getattr(self,field) for field in headings]


  def get_compounds(self, objects=True):
    # Returns either the `CompoundEntry` objects or the `compound` property
    # thereof from the compounds used in this reaction.
    comps = []
    for i in CONFIG.reactant_range():
      comp = getattr(self, "reactant_fk_{}".format(i))
      if comp:
        if objects:
          comps.append(comp)
        else:
          comps.append(comp.compound)
    return comps


  def get_abbrevs(self):
    comps = self.get_compounds()
    return map(lambda comp: comp.abbrev, comps)


  def get_atoms(self, refresh=False):
    if refresh or not self.atoms:
      atoms = {comp.get_atoms(fail_soft=True) for comp in self.get_compounds()}

      # Store `atoms` so that reactions can be efficiently searched by atoms.
      self.atoms = json.dumps(list(atoms))
      self.save()

      return atoms

    else:
      return set(json.loads(self.atoms))


  def refresh(self):
    from DRP.validation import revalidate_datum
    #Store the atoms as a string -- not a set.
    self.get_atoms(refresh=True)

    #Revalidate and save the datum.
    revalidate_datum(self, self.lab_group)






def get_lab_Data(lab_group):
  from DRP.models import get_Lab_Group
  lab_group = get_Lab_Group(lab_group)
  return Data.objects.filter(lab_group=lab_group).order_by("creation_time_dt")


def get_good_rxns(lab_group=None, with_headings=True):
  #Collect a list of all valid data either globally or for a specific lab.

  from DRP.models import get_Lab_Group, convert_QuerySet_to_list

  if lab_group:
    lab_group = get_Lab_Group(lab_group)
    query = get_lab_Data(lab_group).filter(is_valid=True)
  else:
    query = Data.objects.filter(is_valid=True)

  return convert_QuerySet_to_list(query, "Data", with_headings=with_headings)


def get_ref_set(lab_group, reset_cache=True):
  from DRP.cacheFunctions import get_cache, set_cache

  ref_set = get_cache(lab_group, "DATAREFS")
  if not ref_set or reset_cache:
    ref_set = set(get_lab_Data(lab_group).values_list('ref', flat=True))
    set_cache(lab_group, "DATAREFS", ref_set)

  return ref_set


def get_Data_with_compound(comp, lab_group):
  from django.db.models import Q
  from DRP.models import get_compound
  import operator

  if type(comp)==str:
    comp = get_compound(comp, lab_group = lab_group)

  Q_list = [Q(("reactant_fk_{}".format(i),comp)) for i in CONFIG.reactant_range()]
  return Data.objects.filter(reduce(operator.or_, Q_list))


