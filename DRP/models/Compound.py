'''Module containing only the Compound Class'''
from django.db import models, transaction
from MolDescriptor import MolDescriptor
from ChemicalClass import ChemicalClass
from LabGroup import LabGroup
import csv
from chemspipy import ChemSpider
from django.conf import settings
from django.core.exceptions import ValidationError
import importlib

descriptorPlugins = [importlib.import_module(plugin) for plugin in settings.MOL_DESCRIPTOR_PLUGINS]

class CompoundManager(models.Manager):
  '''A custom manager for the Compound Class which permits the creation of entries to and from CSVs'''

  @transaction.atomic
  def fromCsv(self, fileName, labGroup):
    '''Reads a CSV into the database creating objects as a transaction, and returning the resulting queryset of compounds
      (a queryset insures us against very big lists, and allows us to exit the transaction before moving on.
      will get/create for compound classes. This method is all-or-nothing and will fail if one row in the file fails.
      This assumes that the uploaded csv will have headers which map to the names of the fields and that compound classes are
      stored as comma separated lists of the chemicalClass LABEL only.

      Each compound will perform a chemspider-based consistency check on the information it has been created with to ensure
      information is consistent- this throws an ValidationError if it is not.
    '''

    compoundIDsList = []
    cs = ChemSpider(settings.CHEMSPIDER_TOKEN)
    with open(fileName) as f:
      reader = csv.DictReader(f, restkey='restKey')
      rowCount = 0
      errors = []
      for row in reader:
        try:
          rowCount += 1
          if 'chemicalClasses' in row:
            classes = (c.strip() for c in row['chemicalClasses'].split(','))
            chemicalClasses = []
            for c in classes:
              chemicalClass, created = ChemicalClass.objects.get_or_create(label=c)
              chemicalClasses.append(chemicalClass)
          if row.get('CAS') not in ('', None) and row.get('CSID') in ('', None):
            CASResults = cs.simple_search(row['CAS'])
            if len(CASResults) < 1:
              errors.append(ValidationError('CAS Number returned no results from ChemSpider on row %(rowCount)d of uploaded csv.', params={'rowCount':rowCount}))
            elif len(CASResults) == 1:
              row['CSID'] = CASResults[0].csid #a little hacky, but it gets the job done
            else:
              errors.append(ValidationError('CAS number returns more than one ChemSpider ID on row %(rowCount)d of uploaded csv.', params={'rowCount':rowCount}))
          elif row.get('CSID') in ('', None):
            errors.append(ValidationError('No CSID provided on row %(rowCount) of uploaded csv.', params={'rowCount':rowCount}))
          kwargs = {}
          kwargs['CSID'] = row.get('CSID')
          kwargs['abbrev'] = row.get('abbrev')
          kwargs['smiles'] = row.get('smiles')
          kwargs['name'] = row.get('name')
          kwargs['INCHI'] = row.get('INCHI')
          compound = Compound(labGroup = labGroup,  **kwargs)
          compound.full_clean()
          compound.csConsistencyCheck()
          compound.save()
          for chemicalClass in chemicalClasses:
            compound.chemicalClasses.add(chemicalClass)
          compoundIDsList.append(compound.pk)
        except ValidationError as e:
          for message in e.messages:
            errors.append(ValidationError(message + ' on row %(rowCount)d of uploaded csv', params={'rowCount':rowCount}))
      if len(errors) > 0:
        raise ValidationError(errors)
    return compoundIDsList

class Compound(models.Model):
  '''A class for containing data about Compounds used in chemical reactions.
  The assumption is made that all chemicals used are single-species.
  '''
  
  class Meta:
    app_label = "DRP"
    unique_together=(('abbrev', 'labGroup'), ('CSID', 'labGroup'))

  abbrev = models.CharField("Abbreviation", max_length=100)
  '''A local, often nonstandard abbreviation for a compound'''
  name = models.CharField('Name', max_length=300)
  '''Normally the IUPAC name of the compound, however this may not be the most parsable name (which is preferable)'''
  chemicalClasses = models.ManyToManyField(ChemicalClass, verbose_name="Chemical Class")
  '''The class of the compound- examples include Inorganic Salt'''
  CSID = models.PositiveIntegerField('Chemspider ID', null=True)
  '''The chemspider ID for the compound- preferable to the CAS_ID since it is not subject to licensing restrictions'''
  custom = models.BooleanField("Custom", default=False)
  '''This flag denotes whether a compound has been added irrespective of other validation.
  This should be restricted to superusers'''
  INCHI = models.TextField('InCHI key', blank=True, default='')
  '''The Inchi key for a compound- a canonical representation of a molecule which is also unique.'''

  smiles = models.TextField('Smiles', blank=True, default='')
  '''A non-canonical string representation of a molecule which cannot be directly used to test for identity
  but is nevertheless useful for calculating descriptors
  '''

  labGroup = models.ForeignKey(LabGroup, verbose_name="Lab Group")
  '''Tells us whose compound guide this appears in'''

  objects = CompoundManager()

  def csConsistencyCheck(self):
    '''Performs a consistency check of this record against chemspider. Raises a ValidationError on error.'''
    if not self.custom:
      errorList = []
      if self.CSID is None:
        errorList.append('No CSID set')
      cs = ChemSpider(settings.CHEMSPIDER_TOKEN) 
      csCompound = cs.get_compound(self.CSID)
      nameResults = cs.simple_search(self.name)
      if csCompound not in nameResults:
        errorList.append(ValidationError('A compound was consistency checked and was found to have an invalid name', code='invalid_inchi'))
      if self.INCHI == '':
        self.INCHI = csCompound.stdinchi
      elif self.INCHI != csCompound.stdinchi:
        errorList.append(ValidationError('A compound was consistency checked and was found to have an invalid InChi', code='invalid_inchi'))
      if self.smiles == '':
        self.smiles = csCompound.smiles
      elif self.smiles != csCompound.smiles:
        errorList.append(ValidationError('A compound was consistency checked and was found to have an invalid smiles string', code='invalid_smiles'))
      if len(errorList) > 0:
        raise ValidationError(errorList)

  def save(self, commit=True):
    if commit:
      super(Compound, self).save(commit)
      for descriptorPlugin in descriptorPlugins:
        descriptorPlugin.calculate(self) 
    return super(Compound, self).save(commit)
      
