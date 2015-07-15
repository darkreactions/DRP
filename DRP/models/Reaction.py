from django.db import models
from LabGroup import Labgroup
from Descriptor import Descriptor

class Reaction(models.Model):

  class Meta:
    app_label="DRP"

  compounds=models.ManyToManyField(Compound, through="CompoundQuantities")
  temp=models.IntegerField()
  slowCool=models.BooleanField()
  time=models.IntegerField()
  leak=models.BooleanField()
  purity=models.IntegerField()
  notes=models.TextField() 
  labGroup=models.ForeignKey(LabGroup)
  descriptors=models.ManyToManyField(Descriptor, through='DescriptorValue')
