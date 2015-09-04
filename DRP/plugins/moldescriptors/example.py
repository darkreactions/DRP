'''An example molecular descriptor plugin to demonstrate the 'shape' that the API requires.'''
from django.conf import settings
from chemspipy import ChemSpider
from utils import setup
import DRP

_descriptorDict = {
  'mw':{'type':'num','name':'Molecular Weight', 'calculatorSoftware':'chemspider', 'calculatorSoftwareVersion':0, 'maximum':None, 'minimum':0},
  'fs':{'type':'ord', 'name':'Fake size', 'calculatorSoftware':'example.py plugin', 'calculatorSoftwareVersion':0, 'maximum':3, 'minimum':1},
  'N?':{'type':'bool', 'name':'Has Nitrogen', 'calculatorSoftware':'example.py plugin', 'calculatorSoftwareVersion':0},
  'arb':{'type':'cat', 'name':"Phil's arbitrary descriptor", 'calculatorSoftware':'example.py plugin', 'calculatorSoftwareVersion':0, 'permittedValues':('fun', 'dull')}
}
'''A dictionary describing the descriptors available in this module. The key should always be the heading for the descriptor.'''

descriptorDict = setup(_descriptorDict)

def fsValueCalc(mw):
  if mw<50:
    return 1
  elif mw<100:
    return 2
  else:
    return 3

def arbValCalc(compound):
  if compound.pk % 2 == 0:
    return DRP.models.CategoricalDescriptorPermittedValue.objects.get(value='dull', descriptor=descriptorDict['arb'])
  else:
    return DRP.models.CategoricalDescriptorPermittedValue.objects.get(value='fun', descriptor=descriptorDict['arb'])

def calculate(compound):
  '''Calculates the descriptors from this plugin for a compound.
  This should fail silently if a descriptor cannot be calculated for a compound, storing a None value in the
  database as this happens.'''
  cs = ChemSpider(settings.CHEMSPIDER_TOKEN)
  if compound.CSID is None:
    mwValue = DRP.models.NumMolDescriptorValue.objects.get_or_create(compound = compound, descriptor=descriptorDict['mw'], value=None)[0]
    mwValue.save()
    fsValue = DRP.models.OrdMolDescriptorValue.objects.get_or_create(compound=compound, descriptor=descriptorDict['fs'], value=None)[0]
    fsvalue.save()
  else:
    csCompound =  cs.get_compound(compound.CSID)
    mwValue = DRP.models.NumMolDescriptorValue.objects.get_or_create(value=csCompound.molecular_weight, descriptor=descriptorDict['mw'], compound=compound)[0]
    mwValue.save()
    fsValue = DRP.models.OrdMolDescriptorValue.objects.get_or_create(compound=compound, descriptor=descriptorDict['fs'], value=fsValueCalc(mwValue))[0]
    fsValue.save()
  if compound.smiles is None:
    nValue = DRP.models.BoolMolDescriptorValue.objects.get_or_create(compound=compound, descriptor=descriptorDict['N?'], value=None)[0]
  else:
    nValue = DRP.models.BoolMolDescriptorValue.objects.get_or_create(compound=compound, descriptor=descriptorDict['N?'], value=('n' in compound.smiles or 'N' in compound.smiles))[0]
  nValue.save()
  arbValue = DRP.models.CatMolDescriptorValue.objects.get_or_create(compound=compound, descriptor=descriptorDict['arb'], value=arbValCalc(compound))[0]
  arbValue.save()
