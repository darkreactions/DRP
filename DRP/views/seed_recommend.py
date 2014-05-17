# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # # #  Seed Recommendation Views and Functions  # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
from django.http import HttpResponse
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
import django.db

import multiprocessing

from DRP.emailFunctions import email_user
from DRP.retrievalFunctions import *
from DRP.database_construction import *
from DRP.recommendation.seed_rec import constructRecsFromSeed

@login_required
@require_http_methods(["POST"])
def make_seed_recommendations(request):
  u = request.user
  lab_group = u.get_profile().lab_group
  max_seed_calcs = 5

  try:
    pid = request.POST["pid"]
  except:
    return HttpResponse("1") #ERROR: General Failure!

  #Get the active recommendations from the cache.
  active_recs = get_cache(lab_group, "seed_recommendations_active")
  if not active_recs: active_recs = 0

  #Only allow a given number of seed rec calculations at a time.
  if get_cache(lab_group, "seed_recommendations_active")>=max_seed_calcs:
    return HttpResponse("2") #ERROR: Too many active recs at once

  #Start up another seed recommendation calculation if we are able to.
  try:
    #Set the cache and get the appropriate data entry.
    seed = Data.objects.get(id=pid)
    set_cache(lab_group, "seed_recommendations_active", active_recs+1)
   
    #Actually start the new process.
    p = multiprocessing.Process(target=seed_rec_worker, 
		args=(lab_group, seed, u),
		name="{} Seed Rec Worker".format(lab_group.lab_title)
	)
    p.start()

  except Exception as e:
    return HttpResponse("1")
  
  return HttpResponse("0")

def seed_rec_worker(lab_group, seed, user):
  try:
    #Restart the database connection for this new process.
    django.db.close_connection()

    #Actually create new recommendations...
    recList = constructRecsFromSeed(seed.ref) #TODO: As is, this will break if using other Lab Groups.
    #And store them in the database.
    store_new_Recommendation_list(lab_group, recList, seed_source=seed)

    email_body = "The recommendations based on Reaction \"{}\" have finished!".format(seed.ref)
    email_user(user, "Seed Recommendations Ready", email_body)

  except Exception as e:
    print e

    #Email the user that their recommendations failed.
    email_body = "We're very sorry, but the recommendations based on Reaction \"{}\" have failed to be created! Please let us know so that we can fix this!".format(seed.ref)
    email_user(user, "Seed Recommendations Failed!", email_body)

    print "ERROR: Seed recommendation failed! (for \"{}\")".format(lab_group.lab_title)

  #Decrement the number of active recs.
  active_recs = get_cache(lab_group, "seed_recommendations_active")
  
  #In case the cache gets cleared, don't try to subtract from a None type.
  if active_recs:
    set_cache(lab_group, "seed_recommendations_active", active_recs-1)


@login_required
def seed_recommend(request): 
 #Get user data if it exists.
 u = request.user

 fatal_message = ""

 recommendations = get_seed_recs(u.get_profile().lab_group)
 if not recommendations.exists():
   fatal_message = "No recommendations available."   

 return render(request, 'global_page.html', {
  "template":"seed_recommendations",
  "recommendations": recommendations,
  "fatal_message": fatal_message,
 })
