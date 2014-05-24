from retrievalFunctions import *

#Removes any empty Model_Versions that may exist in the database.
def removeEmptyModels():
  print "Beginning empty-model purge."
  models = Model_Version.objects.all()
  i = 0
  for model in models:
    recs_for_model = Recommendation.objects.filter(model_version=model)
    if not recs_for_model.exists():
      model.delete()
      i+=1
  print i
  print "All empty model versions purged ({} total).".format(i)

import datetime
def stringToDateTime(model_name, stringField, dateField):
  #Get the corresponding objects from the database.
  if model_name=="Lab_Member":
    objs = Lab_Member.objects.all()
  elif model_name=="Data":
    objs = Data.objects.all()
  elif model_name=="Recommendation":
    objs = Recommendation.objects.all()
  elif model_name=="Model_Version":
    objs = Model_Version.objects.all()

  i = 0
  skips = 0
  for obj in objs:
    raw_string = getattr(obj, stringField)
    if not raw_string:
      skips += 1
      continue

    dt = datetime.datetime.strptime(raw_string, "%Y-%m-%d %X.%f")
    setattr(obj, dateField, dt)
    obj.save()
    i += 1

  print "{} of {} objects updated ({} skips).".format(i, objs.count(), skips)

