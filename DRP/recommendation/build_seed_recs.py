#!/usr/local/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # Seed Recommendation Gen 'n Store  Worker Process  # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
import sys, os

#Grab the Django settings if they aren't already set.
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

from DRP.settings import BASE_DIR, LOG_DIR
import DRP.models
import django.db

from DRP.emailFunctions import email_user
from DRP.retrievalFunctions import *
from DRP.database_construction import *
from DRP.recommendation.seed_rec import constructRecsFromSeed

#An independent worker process for generating and storing seeds in the database.
def seed_rec_worker(lab_group_id, seed_id, user_id):
  try:
    print "Seed Rec: {} {} {}\n".format(lab_group_id, seed_id, user_id)

    #Restart the database connection for this new process.
    django.db.close_connection()

    #Get the objects from the database (assuming validation has already passed).
    seed = Data.objects.get(id=seed_id)
    lab_group = Lab_Group.objects.get(id=lab_group_id)
    user = User.objects.get(id=user_id)

    #Actually create new recommendations...
    recList = constructRecsFromSeed(seed.ref) #TODO: As is, this will break if using other Lab Groups.
    #And store them in the database.
    store_new_Recommendation_list(lab_group, recList, seed_source=seed)

    email_body = "The recommendations based on Reaction \"{}\" have finished!".format(seed.ref)
    email_user(user, "Seed Recommendations Ready", email_body)

  except Exception as e:
    #Log any errors that might have occurred.
    sys.stderr.write("ERROR: {} {} {}\n".format(lab_group_id, seed_id, user_id))
    sys.stderr.write("{}\n________\n".format(e))
    sys.stderr.flush()

    #Email the user that their recommendations failed.
    email_body = "We're very sorry, but the recommendations based on Reaction \"{}\" could not be created! Please let us know so that we can fix this!".format(seed.ref)
    email_user(user, "Seed Recommendations Failed!", email_body)

  finally: #In the case that emailing fails, always decrement the process count.
    #Decrement the number of active recs.
    active_recs = get_cache(lab_group, "seed_recommendations_active")
  
    #In case the cache gets cleared, don't try to subtract from a None type.
    if active_recs:
      set_cache(lab_group, "seed_recommendations_active", active_recs-1)


if __name__ == "__main__":
  if len(sys.argv) != 4:
    print "You probably want to let the UI handle this..."
    print "python ./this_script.py lab_group_id seed_data_id user_id"
  else:
    seed_rec_worker(sys.argv[1], sys.argv[2], sys.argv[3])
    

