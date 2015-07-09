#!/usr/bin/python

import sys, os

appPath = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.append(appPath)
os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'
