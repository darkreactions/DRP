'''An example molecular descriptor plugin to demonstrate the 'shape' that the API requires.'''
from utils import setup
from django.conf import settings
from chemspipy import ChemSpider
from DRP.models import NumMolDescriptor, BoolMolDescriptor, OrdMolDescriptor, CatMolDescriptor, CatMolDescriptorPermitted  

descriptorDict = {
  'mw':NumMolDescriptor(heading='mw', name='Molecular Weight', calculatorSoftware='chemspider', calculatorSoftwareVersion=0, maximum=None, minimum=0),
  'fs':OrdMolDescriptor(heading='fs', name='Fake size', calculatorSoftware='example.py plugin', calculatorSoftwareVersion=0, maximum=3, minimum=1),
  'N?':BoolMolDescriptor(heading='N?', name='Has Nitrogen', calculatorSoftware='example.py plugin', calculatorSoftwareversion=0)
}
'''A dictionary describing the descriptors available in this module.'''

arb = CatMolDescriptor(heading='arb', name="Phil's arbitrary descriptor", calculatorSoftware='example.py plugin', calculatorSoftwareVersion=0)
arb.save()
fun = CatMolDescriptorPermitted(value='fun', descriptor=arb)
fun.save()
dull =  CatMolDescriptorPermitted(value='dull', descriptor=arb)
dull.save()

descriptorDict['arb'] = arb

for k, v in descriptorDict:
  v.save()

def fsValueCalc(mw):
  if mw<50:
    return 1
  elif mw<100:
    return 2
  else:
    return 3

def calculate(compounds):
  '''Calculates the descriptors from this plugin for each compound in the provided iterable compounds.
  This should fail silently if a descriptor cannot be calculated for a compound, storing a None value in the
  database as this happens.'''
  cs = ChemSpider(settings.CHEMSPIDER_TOKEN)
  for compound in compounds:
    if compound.CSID is None:
      mwValue = NumMolDescriptorValue(compound = compound, descriptor=descriptorDict['mw'], value=None)
      mwValue.save()
      fsValue = OrdMolDescriptor(compound=compound, descriptor=descriptorDict['fs'], value=None)
      fsvalue.save()
    else:
      csCompound =  cs.get_compound(compound.CSID)
      mwValue = csCompound.molecular_weight
      mwValue.save()
      fsValue = OrdMolDescriptor(compound=compound, descriptor=descriptorDict['fs'], value=fsValueCalc(mwValue)
      fsValue.save()
        
