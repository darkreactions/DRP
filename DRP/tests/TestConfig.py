#!/usr/bin/python

import sys, os

path = os.dirname(os.path.realpath(__file__))
if path not in sys.path:
    sys.path.append(path)
