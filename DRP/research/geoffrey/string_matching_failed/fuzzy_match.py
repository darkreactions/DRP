# this was a failure :-(

from munkres import Munkres
from Levenshtein import distance
from sys import argv
import re

# from http://stackoverflow.com/questions/1175208/elegant-python-function-to-convert-camelcase-to-snake-case
first_cap_re = re.compile('(.)([A-Z][a-z]+)')
all_cap_re = re.compile('([a-z0-9])([A-Z])')
def camel_to_snake(name):
    s1 = first_cap_re.sub(r'\1_\2', name)
    return all_cap_re.sub(r'\1_\2', s1).lower()

def strip_and_redorder(s):
    stripped = s.strip()
    snake = camel_to_snake(stripped)
    split_and_sorted = sorted(snake.split('_'))
    sorted_snake = '_'.join(split_and_sorted)
    return sorted_snake
    

# modified from http://stackoverflow.com/questions/26990714/algorithm-for-fuzzy-pairing-of-strings-from-two-lists
if __name__ == '__main__':
    if len(argv) < 3:
        print("Usage: fuzzy_match.py file file")
        print("Finds the best pairing of lines from the two input files")
        print("using the Levenshtein distance and the Hungarian algorithm")
    f1 = argv[1]
    f2 = argv[2]
    
    w1 = [strip_and_redorder(l) for l in open(argv[1]).readlines()]
    w2 = [strip_and_redorder(l) for l in open(argv[2]).readlines()]
    if len(w1) != len(w2):
        if len(w2) > len(w1):
            w1, w2 = w2, w1
        w2.extend([""]*(len(w1)-len(w2)))
    matrix = []
    for i in w1: 
        row = []
        for j in w2:
            d = distance(i.lower(), j.lower())
            if not j:
                d = 10
            row.append(d)
        matrix.append(row)
    m = Munkres()
    max_length = max(len(w) for w in w1)
    for i, j in m.compute(matrix):
        print(("{:<%d}{}" % (max_length+10)).format(w1[i], w2[j]))
