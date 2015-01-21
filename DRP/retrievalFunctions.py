from django.db.models import Q, Max, Min
from models import *

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # DATA  # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Returns a specific datum if it is public or if it belongs to a lab_group.
def get_public_data():
  from DRP.models import Data
  return Data.objects.filter(public=True, is_valid=True).order_by("creation_time_dt")


def get_lab_Data_size(lab_group):
  from DRP.cacheFunctions import get_cache, set_cache
  from DRP.models import get_lab_Data
  size = get_cache(lab_group, "TOTALSIZE")
  if not size:
    size = get_lab_Data(lab_group).count()
    set_cache(lab_group, "TOTALSIZE", size)
  return size


#Get data before/after a specific date (ignoring time).
def filter_by_date(lab_data, raw_date, direction="after"):
  import datetime, dateutil.relativedelta
  # Convert the date input into a usable string. (Date must be given as MM-DD-YY.)
  date = datetime.datetime.strptime(raw_date, "%m-%d-%Y")

  # Get the reactions before/after a specific date.
  if direction.lower() == "after":
    filtered_data = lab_data.filter(creation_time_dt__gte=date)
  else:
    # Add a day to cover any times 00:00-23:59 on a given date.
    date += dateutil.relativedelta.relativedelta(days=1)
    filtered_data = lab_data.filter(creation_time_dt__lte=date)

  return filtered_data


def filter_existing_calcs(data):
  """
  Returns only the data which have calculations.
  """
  from django.db.models import Q

  data = data.filter(~Q(calculations=None))
  data = data.filter(~Q(calculations__contents=""))
  data = data.filter(~Q(calculations__contents="[]"))

  return data


def filter_data(lab_group, query_list):
 from DRP.models import get_lab_Data, get_model_field_names
 from DRP.validation import bool_fields
 import operator

 #Variable Setup
 data = get_lab_Data(lab_group)
 filters = {}
 Q_list = []

  #Collect all the valid search options
 non_reactant_fields = get_model_field_names(unique_only=True)
 foreign_fields = ["user"] #Fields that cannot search by containment.
 reactant_fields = ["reactant","quantity","unit"]
 legal_fields = set(non_reactant_fields+reactant_fields+foreign_fields+["atoms", "public","is_valid"])

 legal_sub_fields = set(get_model_field_names(model="CompoundEntry"))


 #Check the query_list input before performing any database requests.
 for query in query_list:
  try:
   #Make sure values are provided.
   assert query.get(u"field") in legal_fields
   assert query.get(u"match") in {"contain","exact"}
   assert query.get(u"value")

   if query.get(u"sub-field"):
     assert query.get(u"sub-field") in legal_sub_fields

  except:
   raise Exception("One or more inputs is illegal")

 for query in query_list:
  field = query[u"field"]
  sub_field = query[u"sub-field"] if u"sub-field" in query else ""
  match = "__icontains" if query[u"match"] == "contain" else ""
  value = query[u"value"]
  if field in foreign_fields:
   field += "__username" #TODO: Generalize to all fields (make foreign_fields a dict where values are the foreign-field to search).

  #Translate Boolean inputs into Boolean values.
  if field in bool_fields:
   value = True if value.lower()[0] in "1tyc" else False

  #Apply the filter or a Q object with a range of filters.
  if field in reactant_fields:
   or_Qs = []
   for i in CONFIG.reactant_range():
    temp = {field+"_fk_{}".format(i)+"__"+sub_field+match: value}
    or_Qs.append(Q(**temp))

   Q_list.append(reduce(operator.or_, or_Qs))

  elif field=="atoms":
   atom_list = value.split(" ")
   if len(atom_list)>1:
    search_bool = atom_list.pop(-2) #Take the "and" or "or" from the list.
    op = operator.and_ if search_bool == "and" else operator.or_
    atom_Qs = [Q(atoms__contains=atom) for atom in atom_list] #TODO: Test this
    Q_list.append(reduce(op, atom_Qs))
   else:
    Q_list.append(Q(atoms__contains=atom_list[0]))

  else:
   filters[field+match] = value

 #Apply the Q objects and the filters.
 if Q_list:
  data = data.filter(reduce(operator.and_, Q_list))
 if filters:
  data = data.filter(**filters)
 return data


#Either get the valid Data entries for a lab group or get all valid data.
#  Note: Accepts an actual LabGroup object.
def get_valid_data(lab_group=None):
  from models import Data
  if lab_group:
    return Data.objects.filter(is_valid=True, lab_group=lab_group)
  return Data.objects.filter(is_valid=True)



   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # DATACALCS   # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Given some 'data', returns a list of any DataCalc objects that can be gathered.
def expand_data(data, include_lab_info=False, make_new=False, debug=False):
  from DRP.model_building.load_cg import get_cg
  calcList = []

  compound_guide = get_cg()

  for i, datum in enumerate(data):
    if debug and (i%100)==0: print "{}...".format(i)
    try:
      if datum.calculations or make_new:
        calcs = datum.get_calculations_list(include_lab_info=include_lab_info,
                                            preloaded_cg=compound_guide)
        calcList.append(calcs)
    except Exception as e:
      print "(expand_data) {}".format(e)
      pass
  return calcList

#Grabs the "expanded_headers" for the DataCalc objects created in 'parse_rxn.'
def get_expanded_headers():
  from DRP.model_building.rxn_calculator import headers
  return headers


   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # RECOMMENDATIONS # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_recommendations(lab_group):
  from models import Recommendation
  return Recommendation.objects.filter(lab_group=lab_group)

def get_active_recommendations():
  from models import Recommendation
  model = get_latest_ModelStats(active=True)
  return Recommendation.objects.filter(model_version=model)

def get_seed_recs(lab_group, seed_ref=None, show_hidden=False, latest_first=True):
 seed_recs = Recommendation.objects.filter(seeded=True, lab_group=lab_group)

 #If given a seed ref, only yield those recommendations that are seeded from it.
 if seed_ref:
   datum = Data.objects.filter(ref=seed_ref)
   seed_recs = seed_recs.filter(seed=datum)

 if not show_hidden:
   seed_recs = seed_recs.filter(hidden=False)

 if latest_first:
   seed_recs = seed_recs.order_by("-date_dt")

 return seed_recs


def get_recommendations_by_date(lab_group, date = "recent"):
  from DRP.models import Recommendation
  return Recommendation.objects.order_by("-date_dt")


def filter_recommendations(lab_group, query_list):
 from DRP.models import get_model_field_names
 from django.db.models import Q
 import operator

 #Variable Setup
 recs = get_recommendations(lab_group)
 filters = {}
 Q_list = []
  #Collect all the valid search options
 non_reactant_fields = get_model_field_names(model="Recommendation", unique_only=True)
 #Keep track of what field of the ForeignKey should be used to search...
 foreign_fields = {"user":"username",
                   "assigned_user":"username",
                   "seed":"ref"}
 reactant_fields = ["reactant","quantity","unit"]
 legal_fields = set(non_reactant_fields+reactant_fields+["seeded", "saved"]+
                    foreign_fields.keys())

 #Check the query_list input before performing any database requests.
 for query in query_list:
  try:
   #Make sure values are provided.
   assert query.get(u"field") in legal_fields
   assert query.get(u"match") in {"contain","exact"}
   assert query.get(u"value")
  except:
   raise Exception("One or more inputs is illegal")

 for query in query_list:
  field = query[u"field"]
  match = "__icontains" if query[u"match"] == "contain" else ""
  value = query[u"value"]

  #If the field is a ForeignKey, query a field of THAT object.
  if field in foreign_fields:
   field += "__{}".format(foreign_fields[field])

  #Apply the filter or a Q object with a range of filters.
  if field in reactant_fields:
   or_Qs = []
   for i in CONFIG.reactant_range():
    temp = {field+"_{}".format(i)+match: value}
    or_Qs.append(Q(**temp))

   Q_list.append(reduce(operator.or_, or_Qs))
  else:
   filters[field+match] = value

 #Apply the Q objects and the filters.
 if Q_list:
  recs = recs.filter(reduce(operator.and_, Q_list))
 if filters:
  recs = recs.filter(**filters)
 return recs

   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # #  MODEL   # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_usable_models():
  from models import ModelStats
  model_stats = ModelStats.objects.filter(usable=True).order_by("datetime")
  return model_stats


def get_latest_ModelStats(active=True):
  from models import ModelStats
  models = ModelStats.objects.filter(active=active).order_by("-datetime")
  return models.first()


   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # #  LABS AND USERS  # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_lab_users(lab_group):
  lab_members = Lab_Member.filter(lab_group=lab_group)
  users = lab_members.values_list("user", flat=True)
  return users

