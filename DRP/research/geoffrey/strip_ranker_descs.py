#!/usr/bin/env python

import django
from sys import argv
from DRP.models import Descriptor

def strip_desc_headings_threshold(filename, threshold):
    with open(filename, 'rb') as f:
        raw_lines = f.readlines()

    descriptors = []
    found = False
    start_line = "Ranked attributes:"

    for line in raw_lines:
        if found:
            val_num_desc = line.split()
            if len(val_num_desc) == 3:
                val, num, desc = val_num_desc
                if float(val) > threshold:
                    descriptors.append(desc)
        if line.startswith(start_line):
            found = True

    return descriptors

def strip_desc_headings_rank(filename, rank):
    with open(filename, 'rb') as f:
        raw_lines = f.readlines()

    descriptors = []
    found = False
    start_line = "Ranked attributes:"

    for line in raw_lines:
        if found:
            val_num_desc = line.split()
            if len(val_num_desc) == 3:
                val, num, desc = val_num_desc
                descriptors.append(desc)
                if len(descriptors) == rank:
                    break
        if line.startswith(start_line):
            found = True

    return descriptors



def get_descs(headings):
    return [d for d in Descriptor.objects.all() if d.csvHeader in headings]
    # chosen_descs = []
    # for desc in Descriptor.objects.all():
    #     if desc.csvHeader in headings:
    #         chosen_descs.append(desc)
    # return chosen_descs

def get_headings(descs):
    return [d.heading for d in descs]


if __name__ == '__main__':
    django.setup()
    filename = argv[1]
    #rank = int(argv[2])
    threshold = float(argv[2])
    csvHeaders = strip_desc_headings_threshold(filename, threshold)
    #csvHeaders = strip_desc_headings_rank(filename, rank)
    descs = get_descs(csvHeaders)
    headings = get_headings(descs)
    for heading in headings:
        print heading
