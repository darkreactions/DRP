#!/usr/bin/env python
"""
The unit test for the compound Creation form.

These tests assume that presence tests for the form fields work as expected.
"""

import unittest

from django.contrib.auth.models import User

from .baseFormTest import BaseFormTest
from DRP.forms import CompoundForm
from DRP.models import LabGroup, ChemicalClass
from DRP.tests.decorators import createsUser, joinsLabGroup, createsChemicalClass

loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsUser('Aslan', 'old_magic')
@joinsLabGroup('Aslan', 'Narnia')
@createsChemicalClass('Solv', 'Common Solvent')
class NoLabExists(BaseFormTest):
    """Tests that the form doesn't validate when there are no lab groups."""

    def setUpFormData(self):
        """Set up form data."""
        self.formData = {'labGroups': ['5'], 'abbrev': 'etoh', 'name': 'ethanol', 'CAS_ID': '64-17-5', 'CSID': '682',
                         'chemicalClasses': [ChemicalClass.objects.get(label='Solv').pk]}

    def setUp(self):
        """Create a user and a chemical class, then a form."""
        self.user = User.objects.get(username='Aslan')
        self.user.save()
        self.setUpFormData()
        self.form = CompoundForm(self.user, self.formData)


class NoLabForUser(NoLabExists):
    """Tests that the form doesn't validate when there are lab groups but the user is not a member."""

    def setUpFormData(self):
        """Set up form data."""
        super(NoLabForUser, self).setUpFormData()
        self.formData['labGroups'] = [str(self.labGroup.id)]

    def setUp(self):
        """Set up lab group for tests."""
        self.labGroup = LabGroup.objects.makeLabGroup(
            title="LegacyPassTest1", address='1, war drobe, Narnia', email='aslan@example.com', access_code='old_magic')
        self.labGroup.save()
        super(NoLabForUser, self).setUp()

    def tearDown(self):
        """Delete lab group created for this test."""
        self.labGroup.delete()
        super(NoLabForUser, self).tearDown()


class LabForUser(NoLabForUser):
    """Tests that the form validates when the user is a member of the lab. Also checks that the save method works correctly."""

    compound = None

    def setUpFormData(self):
        """Set up form data."""
        super(LabForUser, self).setUpFormData()
        self.labGroup.users.add(self.user)
        self.labGroup.save()

    def test_validation(self):
        """Answer the age old question: Did the test work?."""
        self.validationSucceeds()

    def test_saving(self):
        """Make sure that saving a valid form works."""
        if self.form.is_valid():
            self.compound = self.form.save()
            self.assertIsNotNone(self.compound.id)

    def tearDown(self):
        """Delete objects created for this test."""
        super(LabForUser, self).tearDown()
        if self.compound is not None:
            self.compound.delete()


class LabForBadUser(LabForUser):
    """Tests that the form doesn't validate when the user doesn't choose a lab group."""

    def setUpFormData(self):
        """Get a user and sets its labgroup to an empty string."""
        super(LabForBadUser, self).setUpFormData()
        self.formData['labGroups'] = []

    def test_validation(self):
        """Ensure that the form fails for a user without a lab group."""
        self.validationFails()


class NoCAS(LabForUser):
    """Tests that the form validates with no CAS number provided."""

    def setUpFormData(self):
        """Get data for the form and remove the CAS number."""
        super(NoCAS, self).setUpFormData()
        self.formData['CAS_ID'] = ''


class InconsistentCASName(LabForUser):
    """Test that the form does not validate when the CAS number and name are not consistent."""

    def setUpFormData(self):
        """Create inconsistency between CAS number and name."""
        super(InconsistentCASName, self).setUpFormData()
        self.formData['CAS_ID'] = '290-37-9'
        self.formData['CSID'] = '8904'  # consistent values for pyrazine.

    def test_validation(self):
        """Execute test."""
        self.validationFails()


class InconsistentNameCSID(LabForUser):
    """Tests that the form does not validate when the name and CSID are not consistent."""

    def setUpFormData(self):
        """Create inconsistency between CSID number and name."""
        super(InconsistentNameCSID, self).setUpFormData()
        self.formData['name'] = 'Pyrazine'

    def test_validation(self):
        """Execute test."""
        self.validationFails()


class InconsistentCSIDCAS(LabForUser):
    """Tests that the form does not validate when the CAS number and CSID are not consistent."""

    def setUpFormData(self):
        """Create inconsistency between CAS number and CSID."""
        super(InconsistentCSIDCAS, self).setUpFormData()
        self.formData['CAS_ID'] = '290-37-9'

    def test_validation(self):
        """Execute test."""
        self.validationFails()

suite = unittest.TestSuite([
    loadTests(NoLabExists),
    loadTests(NoLabForUser),
    loadTests(LabForUser),
    loadTests(LabForBadUser),
    loadTests(NoCAS),
    loadTests(InconsistentCASName),
    loadTests(InconsistentNameCSID),
    loadTests(InconsistentCSIDCAS)
])

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
