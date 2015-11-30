'''A module containing decorators which are useful in most test cases for the DRP'''

from DRP.models import Compound, LabGroup, ChemicalClass, License, LicenseAgreement, PerformedReaction, CompoundQuantity, CompoundRole, TestSet, TrainingSet, Descriptor
from DRP.models.rxnDescriptorValues import rxnDescriptorPairs, BoolRxnDescriptorValue, OrdRxnDescriptorValue, NumRxnDescriptorValue, CatRxnDescriptorValue
from DRP.models.rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, CatRxnDescriptor, NumRxnDescriptor
from django.contrib.auth.models import User
from django.conf import settings
from datetime import date, timedelta
import os
import time


def createsRxnDescriptor(heading, descriptorType, options={}):
  '''A class decorator that creates a reaction using pre-existing compounds
     with pre-existing compoundRoles.'''
  def _createsRxnDescriptor(c):
    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    descriptor = None
    for constructor, descriptorVal in rxnDescriptorPairs:
      if descriptorType == constructor.__name__:
        descriptor = constructor()

    if not descriptor:
      error = "Descriptor type \"{}\" unknown to descriptor.".format(descriptorType)
      raise NotImplementedError(error)

    def setUp(self):
      descriptor.heading = heading
      descriptor.name = heading

      for key, val in options.items():
        setattr(descriptor, key, val)

      descriptor.save()

      _oldSetup(self)

    def tearDown(self):
      descriptor.delete()
      _oldTearDown(self)

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _createsRxnDescriptor


def createsPerformedReaction(labTitle, username, compoundAbbrevs, compoundRoles, compoundAmounts, descriptorDict):
  '''A class decorator that creates a reaction using pre-existing compounds
     with pre-existing compoundRoles.'''
  def _createsPerformedReaction(c):
    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    reaction = PerformedReaction()
    compoundQuantities = []
    descriptorVals = []

    def setUp(self):
      labGroup=LabGroup.objects.get(title=labTitle)
      reaction.labGroup = labGroup

      user=User.objects.get(username=username)
      reaction.user = user
      reaction.reference = str(time.time()) #Create a unique reference per reaction.
      reaction.public = False

      reaction.save()

      for abbrev, role, quantity in zip(compoundAbbrevs, compoundRoles, compoundAmounts):
        compound = Compound.objects.get(labGroup=labGroup, abbrev=abbrev)
        compoundRole = CompoundRole.objects.get(label=role)
        compoundQuantity = CompoundQuantity(compound=compound, reaction=reaction,
                                            role=compoundRole, amount=quantity)
        compoundQuantity.save()

        compoundQuantities.append(compoundQuantity)

      for descriptor_heading,val in descriptorDict.items():
        descriptor = Descriptor.objects.filter(heading=descriptor_heading).downcast().next()
        if isinstance(descriptor, BoolRxnDescriptor):
          descriptorVal = BoolRxnDescriptorValue()
        elif isinstance(descriptor, OrdRxnDescriptor):
          descriptorVal = OrdRxnDescriptorValue()
        elif isinstance(descriptor, CatRxnDescriptor):
          descriptorVal = CatRxnDescriptorValue()
        elif isinstance(descriptor, NumRxnDescriptor):
          descriptorVal = NumRxnDescriptorValue()
        else:
          error = "Unknown descriptorValue type for '{}'".format(descriptor)
          raise NotImplementedError(error)

        descriptorVal.descriptor = descriptor
        descriptorVal.value = val
        descriptorVal.reaction = reaction
        descriptorVal.save()

        descriptorVals.append(descriptorVal)

      _oldSetup(self)

    def tearDown(self):
      _oldTearDown(self)

      for cq in compoundQuantities:
        cq.delete()

      for descriptorVal in descriptorVals:
        descriptorVal.delete()

      TrainingSet.objects.filter(reaction=reaction).delete()
      TestSet.objects.filter(reactions__in=[reaction]).delete()
      reaction.delete()


    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _createsPerformedReaction


def createsUser(username, password):
  '''A class decorator that creates a user'''

  def _createsUser(c):

    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    def setUp(self):
      user = User.objects.create_user(username=username, password=password)
      user.save()
      _oldSetup(self)

    def tearDown(self):
      _oldTearDown(self)
      User.objects.filter(username=username).delete()

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _createsUser

def createsCompoundRole(label, description):

  def _createsCompoundRole(c):
    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    role = CompoundRole()

    def setUp(self):
      role.label = label
      role.description = description
      role.save()
      _oldSetup(self)

    def tearDown(self):
      _oldTearDown(self)
      role.delete()

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _createsCompoundRole


def createsCompound(abbrev, csid, classLabel, labTitle, custom=False):

  def _createsCompound(c):

    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    compound = Compound(abbrev=abbrev, CSID=csid, custom=custom)

    def setUp(self):
      compound.labGroup=LabGroup.objects.get(title=labTitle)
      compound.save()
      for c in ChemicalClass.objects.filter(label=classLabel):
        compound.chemicalClasses.add(c)
      compound.save()
      _oldSetup(self)

    def tearDown(self):
      _oldTearDown(self)
      compound.delete()

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _createsCompound

def createsChemicalClass(label, description):
  '''A class decorator that creates a test chemical class for the addition of compounds into the database'''

  def _createsChemicalClass(c):

    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    chemicalClass = ChemicalClass(label=label, description=description)

    def setUp(self):
      chemicalClass.save()
      _oldSetup(self)

    def tearDown(self):
      chemicalClass.delete()
      _oldTearDown(self)

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _createsChemicalClass

def joinsLabGroup(username, labGroupTitle):
  '''A class decorator that creates a test lab group with labGroupTitle as it's title and assigns user identified by
  username to that lab group'''
  def _joinsLabGroup(c):
    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    labGroup = LabGroup(title=labGroupTitle, address='War drobe', email='Aslan@example.com', access_code='new_magic')

    def setUp(self):
      user = User.objects.get(username=username)
      labGroup.save()
      user.labgroup_set.add(labGroup)
      _oldSetup(self)

    def tearDown(self):
      _oldTearDown(self)
      labGroup.delete()

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _joinsLabGroup

def signsExampleLicense(username):
  '''A class decorator that creates a test license and makes the user specified by username sign it on setUp'''
  def _signsExampleLicense(c):

    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    license = License(text='This is an example license used in a test', effectiveDate=date.today() - timedelta(1))

    def setUp(self):
      user = User.objects.get(username=username)
      license.save()
      self.agreement = LicenseAgreement(user=user, text=license)
      self.agreement.save()
      _oldSetup(self)

    def tearDown(self):
      self.agreement.delete()
      license.delete()
      _oldTearDown(self)

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _signsExampleLicense

def loadsCompoundsFromCsv(labGroupTitle, csvFileName):
  '''A class decorator that creates a test set of compounds using the csvFileName, which should be stored in the tests directory resource folder.'''

  def _loadsCompoundsFromCsv(c):

    _oldSetup = c.setUp
    _oldTearDown = c.tearDown

    def setUp(self):
      labGroup = LabGroup.objects.get(title=labGroupTitle)
      compounds = labGroup.compound_set.fromCsv(os.path.join(settings.APP_DIR, 'tests', 'resource', csvFileName))
      for compound in compounds:
        compound.csConsistencyCheck()
        compound.save()
      _oldSetup(self)

    def tearDown(self):
      Compound.objects.all().delete()
      _oldTearDown(self)

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _loadsCompoundsFromCsv
