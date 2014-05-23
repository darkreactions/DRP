from DRP.retrievalFunctions import *

# # # # # # # # # # Seed Recommendations # # # # # # # # # # # # 
def get_seed_rec_worker_list(lab_group):
  active_recs = get_cache(lab_group, "seed_recommendations_active")
  if active_recs: 
    return json.loads(active_recs)
  return []

def cache_seed_rec_worker(lab_group, seed_ref):
  #Either get the cache list or instantiate a new list.
  try:
    active_recs = json.loads(get_cache(lab_group, "seed_recommendations_active"))
  except:
    active_recs = []

  active_recs.append(seed_ref)
  active_rec_json = json.dumps(active_recs)
  set_cache(lab_group, "seed_recommendations_active", active_rec_json)

def remove_seed_rec_worker_from_cache(lab_group, seed):
  #If possible, delete the seed entry from the list of active seed-rec workers.
  active_recs = get_seed_rec_worker_list(lab_group)
  if seed in active_recs: active_recs.remove(seed)

  #Store the new list of active seed-rec workers.
  active_rec_json = json.dumps(active_recs)
  set_cache(lab_group, "seed_recommendations_active", active_rec_json)
