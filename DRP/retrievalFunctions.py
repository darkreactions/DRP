from django.db.models import Q
from models import *

   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # DATA  # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Returns a specific datum if it is public or if it belongs to a lab_group.
def get_public_data():
 #Only show the public data that is_valid.
 return Data.objects.filter(public=True, is_valid=True).order_by("creation_time")

def get_datum_by_ref(lab_group, ref):
 data = get_lab_Data(lab_group)
 datum = data.filter(ref=ref)
 if not datum.exists():
  raise Exception("Datum not found!")
 return datum.first()

def get_lab_Data_size(lab_group):
 size = get_cache(lab_group, "TOTALSIZE")
 if not size:
  size = get_lab_Data(lab_group).count()
  set_cache(lab_group, "TOTALSIZE", size)
 return size

#Get data before/after a specific date (ignoring time).
def get_date_filtered_data(lab_group, raw_date, direction="after", lab_data=None):
 #Convert the date input into a usable string. (Date must be given as MM-DD-YY.)
 date = str(datetime.datetime.strptime(raw_date, "%m-%d-%y"))

 #Only get the data that belongs to a specific lab_group.
 if lab_data:
  lab_data = lab_data.filter(lab_group=lab_group)
 else:
  lab_data = get_lab_Data(lab_group)

 #Get the reactions before/after a specific date.
 if direction.lower() == "after":
  filtered_data = lab_data.filter(creation_time__gte=date_string)
 else:
  filtered_data = lab_data.filter(creation_time__lte=date_string)

 return filtered_data

###Rewrite this to use dict format, Casey. Don't make Paul cry.
def filter_data(lab_group, query_list):
 #Variable Setup:
 lab_data = get_lab_Data(lab_group)
 filters = ""

 #Collect all the valid search options
 non_reactant_fields = get_model_field_names(unique_only=True)
 foreign_fields = ["user"] #Fields that cannot search by containment.
 legal_fields = set(non_reactant_fields+["reactant","quantity","unit","public","is_valid", "atoms"]+foreign_fields)

 #Check the query_list input before performing any database requests.
 try:
  for query in query_list:
   assert query.get(u"field") in legal_fields
   assert query.get(u"match") in {"contain","exact"}
   assert query.get(u"value")
   assert not "\"" in query.get(u"value")

 except:
  raise Exception("One or more inputs is illegal.")
 
 try:
  for query in query_list:
   #Get the query information.
   field = query.get(u"field")
   if field in foreign_fields:
    field = "user__username"
   if query.get(u"match")=="contain" and field not in foreign_fields:
    match = "__icontains"
   else:
    match = ""
   value = query.get(u"value")
 
   if field in list_fields:
    #Check all the reactant/quantity/unit fields.
    Q_obj = ''.join(["Q({}_{}{}=\"{}\")|".format(field, i, match, value) for i in CONFIG.reactant_range()])[:-1]
    filters += ".filter({})".format(Q_obj)
   elif field=="atoms":
    atom_list = value.split(" ")
    if len(atom_list)>1:
     search_bool = atom_list.pop(-2) #Take the "and" or "or" from the list.
     op = "," if search_bool == "and" else "|" #Assign the correct Q operator.
     #Add the atoms  to a Q filter.
     Q_obj = ''.join(["Q(atoms__contains=\"{}\"){}".format(atom, op) for atom in atom_list])[:-1]
    else:
     Q_obj = "Q(atoms__contains=\"{}\")".format(atom_list[0])
    filters += ".filter({})".format(Q_obj)
   else:
    #Translate Boolean inputs into Boolean values.
    if field in bool_fields:
     value = True if value.lower()[0] in "1tyc" else False
     filters += ".filter({}={})".format(field, value)
    else:
     filters += ".filter({}{}=\"{}\")".format(field, match, value)
  data = eval("lab_data"+filters).order_by("creation_time")
  return data

 except Exception as e:
  pass #Security precaution.
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # RECOMMENDATIONS # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_recommendations(lab_group):
 return Recommendation.objects.filter(lab_group=lab_group)

def get_latest_Model_Version(lab_group):
   return Model_Version.objects.filter(lab_group=lab_group, model_type="Recommendation").order_by("-date").first()


def get_recommendations_by_date(lab_group, date = "recent"):
 if date=="recent":
  #Get the most recent version of the model.
  try:
   version = get_latest_Model_Version(lab_group)
   date = version.date
  except Exception as e:
   raise Exception("Could not find any version of the model: {}".format(e))

 #Get the data associated with a specific date.
 try:
  recommendations = get_recommendations(lab_group).filter(date=date).order_by("-score")
 except Exception as e:
  raise Exception("Could not find any version of the model: {}".format(e))

 return recommendations

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
