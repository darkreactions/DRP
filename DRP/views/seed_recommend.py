# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # # #  Seed Recommendation Views and Functions  # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
import django.db
from django.http import HttpResponse
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods

from subprocess import Popen

from DRP.retrievalFunctions import *
from DRP.settings import BASE_DIR, LOG_DIR

@login_required
@require_http_methods(["POST"])
def make_seed_recommendations(request):
  u = request.user
  lab_group = u.get_profile().lab_group
  max_seed_calcs = 5

  try:
    pid = request.POST["pid"]
  except:
    return HttpResponse("1") #1: ERROR: General Failure!

  #Get the active recommendations from the cache.
  active_recs = get_cache(lab_group, "seed_recommendations_active")
  if not active_recs: active_recs = 0

  #Only allow a given number of seed rec calculations at a time.
  if get_cache(lab_group, "seed_recommendations_active")>=max_seed_calcs:
    return HttpResponse("2") #2: ERROR: Too many active recs at once

  #Start up another seed recommendation calculation if we are able to.
  try:
    #Set the cache and get the appropriate data entry.
    set_cache(lab_group, "seed_recommendations_active", active_recs+1)

    #Validate that the recommendation should be started.
    seed = Data.objects.get(id=pid)

    #Get the ids for each object.
    seed_id = seed.id
    lab_id = lab_group.id
    user_id = u.id
   
    #Actually start the new seed-rec construction Process in its own "Pool."
    err_log = open(LOG_DIR+"/seed_recommend/error.log","a")
    act_log = open(LOG_DIR+"/seed_recommend/process.log","a")
    worker_script = BASE_DIR+"/DRP/recommendation/build_seed_recs.py"
    command = "python {} {} {} {}".format(worker_script, lab_id, seed_id, user_id)
    Popen(command.split(), stdout=act_log, stderr=err_log)

  except Exception as e:
    print e
    return HttpResponse("1") #1: ERROR: General Failure!
  
  return HttpResponse("0") #0: Success


@login_required
def seed_recommend(request): 
 #Get user data if it exists.
 u = request.user
 lab_group = u.get_profile().lab_group
 fatal_message = ""

 recommendations = get_seed_recs(lab_group)
 if not recommendations.exists():
   #Get the active recommendations from the cache.
   active_recs = get_cache(lab_group, "seed_recommendations_active")
   if not active_recs: active_recs = 0

   fatal_message = "No recommendations available."
   fatal_message += " (Currently Running: {})".format(active_recs)   

 return render(request, 'global_page.html', {
  "template":"seed_recommendations",
  "recommendations": recommendations,
  "fatal_message": fatal_message,
 })
