from django.db.models import Q, Max, Min
from models import *

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # DATA  # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Returns a specific datum if it is public or if it belongs to a lab_group.
def get_public_data():
 #Only show the public data that is_valid.
 return Data.objects.filter(public=True, is_valid=True).order_by("creation_time_dt")

def get_datum_by_ref(lab_group, ref):
 data = get_lab_Data(lab_group)
 datum = data.filter(ref=ref)
 if not datum.exists():
  raise Exception("Datum not found!")
 return datum[0]

def get_lab_Data_size(lab_group):
 size = get_cache(lab_group, "TOTALSIZE")
 if not size:
  size = get_lab_Data(lab_group).count()
  set_cache(lab_group, "TOTALSIZE", size)
 return size

#Get data before/after a specific date (ignoring time).
def filter_by_date(lab_data, raw_date, direction="after"):
 #Convert the date input into a usable string. (Date must be given as MM-DD-YY.)
 date = datetime.datetime.strptime(raw_date, "%m-%d-%Y")

 #Get the reactions before/after a specific date.
 if direction.lower() == "after":
  filtered_data = lab_data.filter(creation_time_dt__gte=date)
 else:
  filtered_data = lab_data.filter(creation_time_dt__lte=date)

 return filtered_data


def filter_data(lab_group, query_list):
 #Variable Setup
 data = get_lab_Data(lab_group)
 filters = {}
 Q_list = []

  #Collect all the valid search options
 non_reactant_fields = get_model_field_names(unique_only=True)
 foreign_fields = ["user"] #Fields that cannot search by containment.
 reactant_fields = ["reactant","quantity","unit"]
 legal_fields = set(non_reactant_fields+reactant_fields+foreign_fields+["atoms", "public","is_valid"])

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
  if field in foreign_fields:
   field += "__username" #TODO: Generalize to all fields (make foreign_fields a dict where values are the foreign-field to search).

  #Translate Boolean inputs into Boolean values.
  if field in bool_fields:
   value = True if value.lower()[0] in "1tyc" else False

  #Apply the filter or a Q object with a range of filters.
  if field in reactant_fields:
   or_Qs = []
   for i in CONFIG.reactant_range():
    temp = {field+"_{}".format(i)+match: value}
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
  if lab_group:
    return Data.objects.filter(is_valid=True, lab_group=lab_group)
  return Data.objects.filter(is_valid=True)



   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # DATACALCS   # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Given some 'data', returns a list of any DataCalc objects that can be gathered.
def expand_data(data, include_lab_info=False):
  calcList = []
  for datum in data:
    try:
      calcList.append(datum.get_calculations_list(include_lab_info=include_lab_info))
    except Exception as e:
      print e
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
 return Recommendation.objects.filter(lab_group=lab_group)

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

def get_latest_Model_Version(lab_group):
 return Model_Version.objects.filter(lab_group=lab_group, model_type="Recommendation").order_by("-date_dt")[0]


def get_recommendations_by_date(lab_group, date = "recent"):
 if date=="recent":
  #Get the most recent version of the model.
  try:
   version = get_latest_Model_Version(lab_group)
   date = version.date_dt
  except Exception as e:
   raise Exception("Could not find any version of the model: {}".format(e))

 #Get the data associated with a specific date.
 try:
  recommendations = get_recommendations(lab_group).filter(date_dt=date).order_by("-score")
 except Exception as e:
  raise Exception("Could not find any version of the model: {}".format(e))

 return recommendations

def filter_recommendations(lab_group, query_list):
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
   # # # # # # # # # # # # #  COMPOUND GUIDE # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_CG_list(lab_group, headers=True):
 #Variable Setup
 CG = get_lab_CG(lab_group)
 fields =  get_model_field_names(model="CompoundEntry", collect_ignored=True)
 fields.remove("lab_group")
 CG_list = []

 #Apply headers to the list.
 if headers:
  CG_list = fields

 #Get the entries' fields.
 CG_list += [getattr(entry, field) for entry in CG for field in fields]
 return CG_list

def get_CG_list_dict(headers=True):
 labs = Lab_Group.objects.all()
 CG_list_dict = {lab.lab_title:get_CG_list(lab, headers=headers) for lab in labs}
 return CG_list_dict

def get_mass_range(lab_group, abbrev):
 try:
  compounds = get_lab_CG(lab_group).filter(abbrev=abbrev)
  maximum = max([compounds.aggregate(Max("quantity_".format(i))) for i in CONFIG.num_reactants])
  minimum = min([compounds.aggregate(Min("quantity_".format(i))) for i in CONFIG.num_reactants])
  return (minimum, maximuim)
 except:
  raise Exception("Compound or lab not found.")


   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # #  LABS AND USERS  # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_lab_users(lab_group):
  lab_members = Lab_Member.filter(lab_group=lab_group)
  users = lab_members.values_list("user", flat=True)
  return users

