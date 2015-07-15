from django.db import models
from Compound import Compound
from Reaction import Reaction

class CompoundQuantities(models.model):

  class Meta:
    app_label='DRP'

  compound=models.ForeignKey(Compound)
  reaction=models.ForeignKey(Reaction)
  amount=models.FloatField() 
  amountUnit=models.CharField()
