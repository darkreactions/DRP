'''A module containing only the Reaction class'''
from django.db import models
from LabGroup import LabGroup
from Compound import Compound

class Reaction(models.Model):
  '''A base class on which PerformedReactions and RecommendedReactions are built,
  contains common information to each in a table with an automatically
  generated one to one relationship with the subclasses.
  '''

  class Meta:
    app_label="DRP"

  compounds=models.ManyToManyField(Compound, through="CompoundQuantity")
  notes=models.TextField(blank=True) 
  labGroup=models.ForeignKey(LabGroup)
