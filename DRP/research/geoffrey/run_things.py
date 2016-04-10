#!/usr/bin/env python
# import django
# django.setup()
# from DRP.models import *
from sys import argv

fn1 = argv[1]
try:
    fn2 = argv[2]
except IndexError:
    fn2 = fn1

with open(fn1) as f:
    lines = [l.strip() for l in f.readlines()]

lines.sort(key=lambda x: x.lower())

with open(fn2, 'w') as f:
    f.write('\n'.join(lines))
