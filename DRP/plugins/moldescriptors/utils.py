'''A utilities module for helping with molecular descriptor plugins'''
from DRP.models import NumMolDescriptor, BoolMolDescriptor, OrdMolDescriptor, CatMolDescriptor, CatMolDescriptorPermitted  

def setup(descDict):
  returnDict = {}
  for k,v in descDict.items():
    args = v.copy()
    del args['type']
    args['heading']=k
    if v['type'] == 'num':
      returnDict[k] = NumMolDescriptor.objects.get_or_create(**args)
    elif v['type'] == 'bool':
      returnDict[k] = BoolMolDescriptor.objects.get_or_create(**args)
    elif v['type'] == 'ord':
      returnDict[k] = OrdMolDescriptor.objects.get_or_create(**args)
    elif v['type'] == 'cat':
      del args['permittedValues']
      returnDict[k] = CatMolDescriptor(objects.get_or_create(**args)
      returnDict[k].save()
      for permittedValue in v['permittedValues']:
        perm = CatMolDescriptorPermitted.objects.get_or_create(value=permittedValue, descriptor=returnDict[k]) 
        perm.save()
    else:
      raise RuntimeError("Invalid descriptor type provided")
    returnDict[k].save()
