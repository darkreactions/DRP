# Set the Python path so that it has access to the Django settings.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


from DRP.models.Recommendation import *

# CAUTION: Be wary of creating database entries which are not deleted!
# CAUTION: Also be aware of deleting with ForeignKeys! (`clear()` the object!)
# For more: http://stackoverflow.com/questions/8543877/django-models-foreignkey-on-delete-attribute-does-it-really-delete-all-related


def main():
  try:
    #Make a Recommendation
    Recommendation.gather_all_nonsense_recs()
     
  except Exception as e:
    #print "ERROR: {}".format(e)
    error = e 

if __name__=="__main__":
  main()





