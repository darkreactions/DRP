'''Module containing only the Compound Class'''
from django.db import models
from MolDescriptor import MolDescriptor
from ChemicalClass import ChemicalClass
from LabGroup import LabGroup

class Compound(models.Model):
  '''A class for containing data about Compounds used in chemical reactions.
  The assumption is made that all chemicals used are single-species.
  '''
  
  class Meta:
    app_label = "DRP"

  abbrev = models.CharField("Abbreviation", max_length=100)
  '''A local, often nonstandard abbreviation for a compound'''
  name = models.CharField('Name:', max_length=300)
  '''Normally the IUPAC name of the compound, however this may not be the most parsable name (which is preferable)'''
  chemicalClass = models.ManyToManyField(ChemicalClass, verbose_name="Chemical Class")
  '''The class of the compound- examples include Inorganic Salt'''
  CSID = models.PositiveIntegerField('Chemspider ID', null=True)
  '''The chemspider ID for the compound- preferable to the CAS_ID since it is not subject to licensing restrictions'''
  custom = models.BooleanField("Custom", default=False)
  '''This flag denotes whether a compound has been added irrespective of other validation.
  This should be restricted to superusers'''
  INCHI = models.TextField('InCHI key', blank=True, default='')
  '''The Inchi key for a compound- a canonical representation of a molecule which is also unique.'''

  smiles= models.TextField('Smiles', blank=True, default='')
  '''A non-canonical string representation of a molecule which cannot be directly used to test for identity
  but is nevertheless useful for calculating descriptors
  '''

  descriptors = models.ManyToManyField(MolDescriptor, through='MolDescriptorValue')
  '''A link to descriptors which have been calculated for this compound. Values for the descriptors are found
  on the MolDescriptorValue model.
  '''

  labGroup = models.ForeignKey(LabGroup)
  '''Tells us whose compound guide this appears in'''
