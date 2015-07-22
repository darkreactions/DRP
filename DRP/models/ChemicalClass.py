'''A module containing only the ChemicalClass class'''
from django.db import models
from DRP.settings import CHEMICAL_CLASS_CHOICES

class ChemicalClass(models.Model):

  class Meta:
    app_label='DRP'

  label=models.CharField(max_length=30, unique=True, error_messages={'unique':'A chemical class with this label already exists'}, choices=CHEMICAL_CLASS_CHOICES)
  description=models.CharField(max_length=20)
