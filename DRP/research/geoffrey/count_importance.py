#!/usr/bin/env python
import subprocess
import sys
import os
import tempfile
import json

fd, name = tempfile.mkstemp()
os.close(fd)
os.system("./parse_j48.py %s > %s" % (sys.argv[1], name))
j = json.load(open(name))
os.unlink(name)

importance = {}


def walk_tree(j, feature_set=set()):
    feature, children = j
    for child in children:
        if len(child) == 3:
            walk_tree(child[2], feature_set.union(set([feature])))
        elif len(child) == 4:
            for leaf_feature in feature_set:
                importance[leaf_feature] = importance.get(
                    leaf_feature, 0) + child[3][0] + child[3][1]

walk_tree(j)

rank = sorted(importance.items(), key=lambda k_v: k_v[1], reverse=True)
mx = float(rank[0][1])

max_width = max(map(len, importance.keys()))
print "{feature:{width}} {prob:10} {rank:>5}".format(feature="Feature", width=max_width + 5, prob="Importance", rank="Rank")

for r, (feature, count) in enumerate(rank):
    # print "%s,\t %0.3f,\t %s" % (feature, count / float(mx), r+1)
    print "{feature:<{width}} {prob:10.3f} {rank:5}".format(feature=feature, width=max_width + 5, prob=count / mx, rank=r + 1)
