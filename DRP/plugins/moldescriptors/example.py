'''An example molecular descriptor plugin to demonstrate the 'shape' that the API requires.'''
#from django.conf import settings
from chemspipy import ChemSpider
from DRP.models import NumMolDescriptorValue, BoolMolDescriptorValue, OrdMolDescriptorValue, CatMolDescriptorValue, CatMolDescriptorPermitted  
from utils import setup

_descriptorDict = {
  'mw':{'type':'num','name':'Molecular Weight', 'calculatorSoftware':'chemspider', 'calculatorSoftwareVersion':0, 'maximum':None, 'minimum':0},
  'fs':{'type':'ord', 'name':'Fake size', 'calculatorSoftware':'example.py plugin', 'calculatorSoftwareVersion':0, 'maximum':3, 'minimum':1),
  'N?':{'type':'bool', 'name':'Has Nitrogen', 'calculatorSoftware':'example.py plugin', 'calculatorSoftwareversion':0}
  'arb':{'type':'cat', 'name':"Phil's arbitrary descriptor", 'calculatorSoftware':'example.py plugin', 'calculatorSoftwareVersion':0, 'permittedValues':('fun', 'dull')}
}
'''A dictionary describing the descriptors available in this module. The key should always be the heading for the descriptor.'''

descriptorDict['arb'] = arb

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
    return CatMolDescriptorPermitted.get(value='dull', descriptor=descriptorDict['arb'])
  else:
    return CatMolDescriptorPermitted.get(value='fun', descriptor=descriptorDict['arb'])

def calculate(compound):
  '''Calculates the descriptors from this plugin for a compound.
  This should fail silently if a descriptor cannot be calculated for a compound, storing a None value in the
  database as this happens.'''
#  cs = ChemSpider(settings.CHEMSPIDER_TOKEN)
  if compound.CSID is None:
    mwValue = NumMolDescriptorValue(compound = compound, descriptor=descriptorDict['mw'], value=None)
    mwValue.save()
    fsValue = OrdMolDescriptorValue(compound=compound, descriptor=descriptorDict['fs'], value=None)
    fsvalue.save()
  else:
    csCompound =  cs.get_compound(compound.CSID)
    mwValue = csCompound.molecular_weight
    mwValue.save()
    fsValue = OrdMolDescriptor(compound=compound, descriptor=descriptorDict['fs'], value=fsValueCalc(mwValue))
    fsValue.save()
  if compound.smiles is None:
    nValue = BoolMolDescriptorValue(compound=compound, descriptor=descriptorDict['N?'], value=None)
  else:
    nValue = BoolMolDescriptorValue(compound=compound, descriptor=descriptorDict['N?'], value=('n' in smiles or 'N' in smiles))
  nValue.save()
  arbValue = CatMolDescriptorValue(compound=compound, descriptor=descriptorDict['arb'], value=arbValCalc(compound))
  arbValue.save()
