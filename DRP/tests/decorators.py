'''A module containing decorators which are useful in most test cases for the DRP'''

from chemspipy import ChemSpider
from DRP.models import Compound, LabGroup, ChemicalClass, License, LicenseAgreement
from django.contrib.auth.models import User
from django.conf import settings
from datetime import date, timedelta
import os

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
      User.objects.filter(username=username).delete()
      _oldTearDown(self)

    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _createsUser

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
  '''A class decorators that creates a test set of compounds using the csvFileName, which should be stored in the tests directory resource folder.'''
  
  def _loadsCompoundsFromCsv(c):

    _oldSetup = c.setUp
    _oldTearDown = c.tearDown
    labGroup = LabGroup.objects.get(title=labGroupTitle)

    def setUp(self):
      compounds = labGroup.compounds.fromCsv(os.path.join(settings.APP_PATH, 'tests', 'resource', csvFileName)
      for compound in compounds:
        compound.save()

    def tearDown(self):
      Compound.objects.all().delete()
    
    c.setUp = setUp
    c.tearDown = tearDown
    return c
  return _loadsCompoundsFromCsv
