#!/usr/bin/env python
from sys import argv

filenames = argv[1:]

h_set = set()

for fn in filenames:
    with open(fn) as f:
        new_h = set([l.strip() for l in f.readlines() if l.strip()])
        h_set = h_set.union(new_h)

print '\n'.join(sorted(h_set, key=lambda x: x.lower())),
