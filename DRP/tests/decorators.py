'''A module containing decorators which are useful in most test cases for the DRP'''

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
