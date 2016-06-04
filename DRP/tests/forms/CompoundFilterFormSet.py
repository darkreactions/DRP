#!/usr/bin/env python
'''Test classes for the compound filter form'''

import django
django.setup()

import unittest
from DRP.forms import FilterForm, CompoundFilterFormSet
from BaseFormTest import BaseFormTest
from DRP.models import LabGroup, ChemicalClass
from django.conf import settings
from django.contrib.auth.models import User
from DRP.tests.decorators import createsCompound, createsUser, joinsLabGroup, createsChemicalClass, loadsCompoundsFromCsv
from DRP.tests.DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@createsChemicalClass("Amine", "descr")
@createsChemicalClass("Inorg", "descr")
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
class CompoundFilterSetSucceed(BaseFormTest):
    '''Three different forms searching for different compounds, should all succeed'''

    def setUpFormData(self):
        self.formData = {}
        self.formData["form-TOTAL_FORMS"] = 3
        self.formData["form-INITIAL_FORMS"] = 0
        self.formData["form-MAX_NUM_FORMS"] = 10000
        # first form
        self.formData["form-0-abbrevs"] = "2-amep"
        self.formData["form-0-name"] = "2-aminomethyl-1-ethylpyrrolidine"
        self.formData["form-0-chemicalClasses"] = [self.chemicalClass.pk]
        self.formData["form-0-CSID"] = "104820"
        self.formData["form-0-INCHI"] = "InChI=1S/C7H16N2/c1-2-9-5-3-4-7(9)6-8/h7H,2-6,8H2,1H3"
        self.formData["form-0-smiles"] = "CCN1CCCC1CN"
        self.formData["form-0-labGroup"] = self.labgroup.pk
        self.formData["form-0-js_active"] = False
        # second form
        self.formData["form-1-abbrevs"] = "hmta"
        self.formData["form-1-name"] = "1,3,5,7-tetraazaadamantane"
        self.formData["form-1-chemicalClasses"] = [self.chemicalClass.pk]
        self.formData["form-1-CSID"] = "3959"
        self.formData["form-1-INCHI"] = "InChI=1S/C6H12N4/c1-7-2-9-4-8(1)5-10(3-7)6-9/h1-6H2"
        self.formData["form-1-smiles"] = "C1N2CN3CN1CN(C2)C3"
        self.formData["form-1-labGroup"] = self.labgroup.pk
        self.formData["form-1-js_active"] = False
        # third form
        self.formData["form-2-abbrevs"] = "sov"
        self.formData["form-2-name"] = "Sodium Orthovanadate"
        self.formData["form-2-chemicalClasses"] = [self.chemicalClass2.pk]
        self.formData["form-2-CSID"] = "55575"
        self.formData["form-2-INCHI"] = "InChI=1S/3Na.4O.V/q3*+1;;3*-1;"
        self.formData["form-2-smiles"] = "[O-][V](=O)([O-])[O-].[Na+].[Na+].[Na+]"
        self.formData["form-2-labGroup"] = self.labgroup.pk
        self.formData["form-2-js_active"] = False

    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora')
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.chemicalClass = ChemicalClass.objects.get(label='Amine')
        self.chemicalClass2 = ChemicalClass.objects.get(label='Inorg')
        self.setUpFormData()
        self.form = CompoundFilterFormSet(user=self.user, labGroup=self.labgroup, data=self.formData)

    def test_validation(self):
        '''test that this test returned correct compounds'''
        self.validationSucceeds()

    def test_is_empty(self):
        self.assertTrue(self.form.is_valid(), str(self.form.errors))
        self.assertFalse(self.form.is_empty())

    def test_cleaned_data(self):
        self.assertTrue(self.form.is_valid())
        cleaned = self.form.cleaned_data
        self.assertEqual(len(cleaned), 3, "Got more than one result from the filter form, or cleaned more than the results returned from the form")

    def test_fetch(self):
        self.assertTrue(self.form.is_valid())
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 3, "{} results found.".format(fetched.count()))
        self.assertEqual(fetched.get(abbrev="2-amep").abbrev, "2-amep", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))
        self.assertEqual(fetched.get(abbrev='hmta').abbrev, "hmta", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[1].abbrev))
        self.assertEqual(fetched.get(abbrev="sov").abbrev, "sov", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[2].abbrev))


suite = unittest.TestSuite([
    loadTests(CompoundFilterSetSucceed),
])

if __name__ == '__main__':
    runTests(suite)
