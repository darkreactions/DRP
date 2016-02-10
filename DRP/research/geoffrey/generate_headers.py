#!/usr/bin/env python

from DRP.models import Descriptor

def print_headers():
  for descriptor in Descriptor.objects.all():
    print descriptor.heading

if __name__=='__main__':
  print_headers()
