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
