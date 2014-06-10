# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # # #  Seed Recommendation Views and Functions  # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
import django.db
from django.http import HttpResponse
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods

import json
from subprocess import Popen

from DRP.retrievalFunctions import *
from DRP.settings import BASE_DIR, LOG_DIR
from DRP.cacheFunctions import *

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
  active_workers = get_seed_rec_worker_list(lab_group)

  #Only allow a given number of seed rec calculations at a time.
  if len(active_workers)>=max_seed_calcs:
    return HttpResponse("2") #2: ERROR: Too many active recs at once

  #Start up another seed recommendation calculation if we are able to.
  try:
    #Validate that the recommendation should be started.
    seed = Data.objects.get(id=pid)

    #Get the ids for each object.
    seed_id = seed.id
    lab_id = lab_group.id
    user_id = u.id

    #Set the cache and get the appropriate data entry.
    cache_seed_rec_worker(lab_group, seed.ref)
   
    #Actually start the new seed-rec construction process to build recs.
    err_log = open(LOG_DIR+"/seed_recommend/error.log","a")
    act_log = open(LOG_DIR+"/seed_recommend/process.log","a")
    worker_script = BASE_DIR+"/DRP/recommendation/build_seed_recs.py"
    command = "python {} {} {} {}".format(worker_script, lab_id, seed_id, user_id)
    #Log to the files above and make the worker independent of the parent process.
    Popen(command.split(), stdout=act_log, stderr=err_log, close_fds=True)

  except Exception as e:
    print e
    return HttpResponse("1") #1: ERROR: General Failure!
  
  return HttpResponse("0") #0: Success


@login_required
def check_seed_worker_cache(request):
  #Get any seed_workers that are in the cache.
  lab_group = request.user.get_profile().lab_group
  active_workers = get_seed_rec_worker_list(lab_group)

  #And send back the updated UI.
  return render(request, 'seed_recommendations_cache_display.html', {
    "currently_running": active_workers,
   })


@login_required
def seed_recommend(request): 
   #Get user data if it exists.
   u = request.user
   lab_group = u.get_profile().lab_group
   fatal_message = ""
  
   try:
     #Either get all of the recommendations or just the non-hidden ones.
     show_hidden = request.GET.get("show_hidden")=="True"
     recommendations = get_seed_recs(lab_group, show_hidden=show_hidden)[:200]
     #TODO:Ideally, pagify this, Future Casey.
  
     #Get the active recommendations from the cache.
     active_workers = get_seed_rec_worker_list(lab_group)
  
     assert recommendations.exists()
   except Exception as e:
     fatal_message = "No recommendations available."
     active_workers = []
     recommendations = []
  
   return render(request, 'global_page.html', {
    "template":"seed_recommendations",
    "recommendations": recommendations,
    "fatal_message": fatal_message,
    "currently_running": active_workers,
   })
  
