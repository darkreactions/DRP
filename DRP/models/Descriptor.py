'''A module containing only the Descriptor class'''
from django.db import models

class Descriptor(models.Model):
  '''A class which describes a descriptor- a value which describes a system such as a compound or a reaction'''
  
  class Meta:
    app_label='DRP'

  heading=models.CharField(max_length=200, unique=True, error_messages={'unique':'This descriptor is already registered, or another descriptor already has this title.'})
  '''A short label which is given to a description. No constraints currently exist, but this may be tweaked later to
  enforce MS-excel style CSV compatibility
  '''
  name=models.CharField('Full name', max_length=300)
  kind=models.CharField('Kind', max_length=20, choices=(('Cat', 'Categorical'), ('Ord', 'Ordinal'), ('Num', 'Numerical'), ('Bool', 'Boolean')))
  '''The kind of descriptor changes the type of data which is applicable for a descriptor value (these are stores in the DescriptorValue relationship class'''
