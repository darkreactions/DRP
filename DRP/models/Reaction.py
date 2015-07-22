'''A module containing only the Reaction class'''
from django.db import models
from LabGroup import Labgroup
from Descriptor import Descriptor

class Reaction(models.Model):
  '''A base class on which PerformedReactions and RecommendedReactions are built,
  contains common information to each in a table with an automatically
  generated one to one relationship with the subclasses.
  '''

  class Meta:
    app_label="DRP"

  compounds=models.ManyToManyField(Compound, through="CompoundQuantity")
  temp=models.IntegerField()
  slowCool=models.BooleanField()
  time=models.IntegerField()
  leak=models.BooleanField()
  purity=models.IntegerField()
  notes=models.TextField() 
  labGroup=models.ForeignKey(LabGroup)
  descriptors=models.ManyToManyField(Descriptor, through='DescriptorValue')
