'''A module containing Classes permitting the representation of molecular descriptors'''
from django.db import models

class MolDescriptor(models.Model):
  '''An abstract class which describes a descriptor- a value which describes a system such as a compound or a reaction'''
  
  class Meta:
    app_label='DRP'
    verbose_name = 'Molecular Descriptor'
    unique_together = ('heading','calculatorSoftware','calculatorSoftwareVersion')

  heading=models.CharField(max_length=200, unique=True, error_messages={'unique':'This descriptor is already registered, or another descriptor already has this title.'})
  '''A short label which is given to a description. No constraints currently exist, but this may be tweaked later to
  enforce MS-excel style CSV compatibility
  '''
  name=models.CharField('Full name', max_length=300)
  kind=models.CharField('Kind', max_length=20, choices=(('Cat', 'Categorical'), ('Ord', 'Ordinal'), ('Num', 'Numerical'), ('Bool', 'Boolean')))
  '''The kind of descriptor changes the type of data which is applicable for a descriptor value (these are stores in the DescriptorValue relationship class'''
  calculatorSoftware=models.CharField(max_length=100)
  calculatorSoftwareVersion=models.CharField(max_length=20)

class CatMolDescriptor(MolDescriptor):
  '''A class which describes a categorical molecular descriptors'''

  class Meta:
    app_label='DRP'
    verbose_name= 'Categorical Molecular Descriptor'

class OrdMolDescriptor(MolDescriptor):
  '''A class which represents an ordinal descriptor'''
  
  class Meta:
    verbose_name= 'Ordinal Molecular Descriptor'
    app_label='DRP'

  maximum=models.IntegerField()
  minimum=models.IntegerField()

class NumMolDescriptor:
  '''A class which represents a numerical descriptor'''

  class Meta:
    app_label='DRP'
    verbose_name= 'Numerical Molecular Descriptor'

  maximum=models.FloatField()
  minimum=models.FloatField()

class BoolMolDescriptor(MolDescriptor):
  '''A class which represents a boolean descriptors'''

  class Meta:
    app_label='DRP'
    verbose_name= 'Boolean Molecular Descriptor'

class CatMolDescriptorPermitted(models.Model):
  '''A class which represents the permitted values for a categorical descriptor'''

  class Meta:
    app_label = "DRP"
    verbose_name= 'Permitted Categorical Descriptor Value'

  descriptor=models.ForeignKey(CatMolDescriptor, related_name='permittedValues')
