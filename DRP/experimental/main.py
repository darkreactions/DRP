#!/usr/bin/env python

#Grab the Django settings if they aren't already set.
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
  sys.path = [django_path] + sys.path
  os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

def run():
  from DRP.experimental.model_building import build_previous_model
  date = "06-01-2014"
  title = "Retrogenerated {}".format(date)
  description = "A model retrogenerated from the data available on {}".format(date)
  build_previous_model(title, description, date)

if __name__=="__main__":
  run()
