from os import path
import sys

drpTestPath = path.dirname(os.realpath(__FILE__))
drpPath = path.dirname(drpTestPath)

if drpPath not in sys.path:
    sys.path.append(drpPath)
