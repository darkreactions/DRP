
'''Module containing only the Compound Class'''
from django.db import models, transaction
from MolDescriptor import MolDescriptor, BoolMolDescriptor, NumMolDescriptor, CatMolDescriptor, OrdMolDescriptor
from MolDescriptorValue import BoolMolDescriptorValue, NumMolDescriptorValue, CatMolDescriptorValue, OrdMolDescriptorValue
from ChemicalClass import ChemicalClass
from LabGroup import LabGroup
import csv
from querysets import CsvModel, CsvQuerySet, ArffQuerySet
from chemspipy import ChemSpider
from django.conf import settings
from django.core.exceptions import ValidationError
from itertools import chain
import importlib
from collections import OrderedDict

descriptorPlugins = [importlib.import_module(plugin) for plugin in settings.MOL_DESCRIPTOR_PLUGINS] #this prevents a cyclic dependency problem

class CompoundQuerySet(CsvQuerySet, ArffQuerySet):
  '''A specialised queryset for outputting Compounds in specific formats'''

  def __init__(self, **kwargs):
    '''Initialises teh queryset'''
    kwargs.pop('model', None)
    super(CompoundQuerySet, self).__init__(Compound, **kwargs)

  def maxChemicalClassCount(self):
    '''Gives a count of the maximum number of chemical classes associated with a compound in this queryset'''
    m = self.annotate(chemicalClassCount=models.Count('chemicalClasses')).aggregate(max=models.Max('chemicalClassCount'))['max']
    return 0 if m is None else m
    
  @property
  def csvHeaders(self):
    '''Generates the header row information for the CSV'''
    headers = super(CompoundQuerySet, self).csvHeaders
    m = Compound.objects.all().maxChemicalClassCount()
    headers += ['chemicalClass_{}'.format(x+1) for x in range(0, m)]
    return headers

  @property
  def arffHeaders(self):
    '''generates headers for the arff file'''
    headers = super(CompoundQuerySet, self).arffHeaders
    m = Compound.objects.all().maxChemicalClassCount()
    for x in range(0, m):
      label = 'chemicalClass_{0}'.format(x+1)
      headers[label] = '@attribute {} {{{}}}'.format(label, ','.join(str(chemicalClass) for chemicalClass in ChemicalClass.objects.all()))
    return headers

  @property
  def expandedArffHeaders(self):
    headers = self.arffHeaders
    headers.update(OrderedDict(((d.csvHeader, d.arffHeader) for d in self.descriptors())))
    return headers

  @property
  def expandedCsvHeaders(self):
    '''Generates the expanded header for the csv'''
    return self.csvHeaders + [d.csvHeader for d in self.descriptors()]

  def descriptors(self):
    '''returns the descriptor which have relationship to the queryset'''
    return chain(
        BoolMolDescriptor.objects.filter(
          boolmoldescriptorvalue__in=BoolMolDescriptorValue.objects.filter(
            compound__in=self
          )
        ),
        NumMolDescriptor.objects.filter(
          nummoldescriptorvalue__in=NumMolDescriptorValue.objects.filter(
            compound__in=self
          )
        ),
        OrdMolDescriptor.objects.filter(
          ordmoldescriptorvalue__in=OrdMolDescriptorValue.objects.filter(
            compound__in=self
          )
        ),
        CatMolDescriptor.objects.filter(
          catmoldescriptorvalue__in=CatMolDescriptorValue.objects.filter(
            compound__in=self
          )
        )
      )
    

class CompoundManager(models.Manager):
  '''A custom manager for the Compound Class which permits the creation of entries to and from CSVs'''

  use_for_related_fields = True #NOTE:This doesn't actually work, but no-one's sure which way django is going to jump on this.

  def get_queryset(self):
    return CompoundQuerySet()

  def fromCsv(self, fileName, labGroup=None):
    '''Reads a CSV into the creating objects, returning a list of compounds which have not yet been saved.
      This assumes that the uploaded csv will have headers which map to the names of the fields and that compound classes are
      stored as comma separated lists of the chemicalClass LABEL only.

      Each compound will perform a chemspider-based consistency check on the information it has been created with to ensure
      information is consistent- this throws an ValidationError if it is not.
    '''
    if labGroup is None and hasattr(self, 'instance'):
      #we presume that if this is being called without a labgroup that's because this manager belongs to a lab group
      labGroup = self.instance

    compoundsList = []
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
            errors.append(ValidationError('No CSID provided on row %(rowCount)d of uploaded csv.', params={'rowCount':rowCount}))
          kwargs = {}
          kwargs['CSID'] = row.get('CSID')
          kwargs['abbrev'] = row.get('abbrev')
          kwargs['smiles'] = row.get('smiles')
          kwargs['name'] = row.get('name')
          kwargs['INCHI'] = row.get('INCHI')
          compound = Compound(labGroup = labGroup,  **kwargs)
          for chemicalClass in chemicalClasses:
            compound.lazyChemicalClasses.append(chemicalClass)
          compoundsList.append(compound)
        except ValidationError as e:
          for message in e.messages:
            errors.append(ValidationError(message + ' on row %(rowCount)d of uploaded csv', params={'rowCount':rowCount}))
      if len(errors) > 0:
        raise ValidationError(errors)
    return compoundsList

class Compound(CsvModel):
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

  def __init__(self, *args, **kwargs):
    super(Compound, self).__init__(*args, **kwargs)
    self.lazyChemicalClasses= []

  def csConsistencyCheck(self):
    '''Performs a consistency check of this record against chemspider. Raises a ValidationError on error.'''
    if not self.custom:
      errorList = []
      cs = ChemSpider(settings.CHEMSPIDER_TOKEN) 
      if self.CSID is None or self.CSID is '':
        raise ValidationError('No CSID set', 'no_csid')
      else:
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

  def save(self, *args, **kwargs):
    super(Compound, self).save(*args, **kwargs)
    for lcc in self.lazyChemicalClasses: #coping mechanism for compounds loaded from csv files; not to be used by other means
      self.chemicalClasses.add(lcc)
    for descriptorPlugin in descriptorPlugins:
      descriptorPlugin.calculate(self) 
    
  @property  
  def descriptorValues(self):
    '''Returns a union of the descriptor values that apply to this object. Hijacks the method from the CompoundQuerySet class'''
    return set(chain(self.catmoldescriptorvalue_set.all(), self.nummoldescriptorvalue_set.all(), self.ordmoldescriptorvalue_set.all(), self.boolmoldescriptorvalue_set.all()))

  @property
  def values(self):
    '''outputs a dict of values suitable for use by a csv.DictWriter'''
    d =  super(Compound, self).values 
    i = 1
    for c in self.chemicalClasses.all():
      d['chemicalClass_{}'.format(i)] = c
      i+=1
    return d

  @property
  def expandedValues(self):
    '''outputs a dict of values suitable for use by a csv.DictWriter - includes molecular descriptors'''
    res = self.values.copy()
    res.update({descriptorValue.descriptor.csvHeader:descriptorValue.value for descriptorValue in self.descriptorValues})
    return res 
