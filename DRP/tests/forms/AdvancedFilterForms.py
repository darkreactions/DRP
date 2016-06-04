#!/usr/bin/env python
'''Test classes for the compound filter form'''

import unittest 
from DRP.forms import FilterForm
from DRP.forms.compound.filterforms import NumericFilterForm, CategoryFilterForm, BooleanFilterForm, OrdinalFilterForm
from DRP.models import CatMolDescriptor, BoolMolDescriptor, NumMolDescriptor, OrdMolDescriptor
from BaseFormTest import BaseFormTest
from DRP.models import LabGroup, ChemicalClass
from django.conf import settings
from django.contrib.auth.models import User
from DRP.tests.decorators import createsCompound, createsUser, joinsLabGroup, createsChemicalClass, loadsCompoundsFromCsv
from DRP.tests.DRPTestCase import DRPTestCase, runTests
from DRP.models import CategoricalDescriptorPermittedValue
loadTests = unittest.TestLoader().loadTestsFromTestCase


@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
class NumericFilterFormGT(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = NumMolDescriptor.objects.get(heading="length").pk
        self.formData["operator"] = 'gt' 
        self.formData["value"] = 12

    def setUp(self): 
        '''Creates a user, then a form'''
        self.setUpFormData()
        self.form = NumericFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct compounds'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid() 
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["descriptor"], NumMolDescriptor.objects.get(heading="length"), "Got a different result, {}, for the descriptor than expected".format(cleaned["descriptor"]))
        self.assertEqual(cleaned["operator"], "gt", "Got a different result, {}, for the operator than expected".format(cleaned["descriptor"]))
        self.assertEqual(float(cleaned["value"]), 12, "Got a different result, {}, for the value than expected".format(cleaned["value"]))

    def test_is_empty(self):
        self.form.is_valid() 
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got {} objects ( {} ) instead of 1".format(fetched.count(), [compound for compound in fetched]))
        self.assertEqual(fetched[0].compound.abbrev, "hmta", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].compound.abbrev))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
class NumericFilterFormGE(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = NumMolDescriptor.objects.get(heading="length").pk
        self.formData["operator"] = 'ge' 
        self.formData["value"] = 18

    def setUp(self): 
        '''Creates a user, then a form'''
        self.setUpFormData()
        self.form = NumericFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct compounds'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["descriptor"], NumMolDescriptor.objects.get(heading="length"), "Got a different result, {}, for the descriptor than expected".format(cleaned["descriptor"]))
        self.assertEqual(cleaned["operator"], "ge", "Got a different result, {}, for the operator than expected".format(cleaned["operator"]))
        self.assertEqual(float(cleaned["value"]), 18, "Got a different result,{}, for the value than expected".format(cleaned["value"]))

    def test_is_empty(self):
        self.form.is_valid() 
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 2, "Got {} results intead of 2".format(fetched.count()))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
class OrdinalFilterFormLT(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = OrdMolDescriptor.objects.get(heading="fs").pk
        self.formData["operator"] = "lt"
        self.formData["value"] = 3

    def setUp(self): 
        self.setUpFormData()
        self.form = OrdinalFilterForm(data=self.formData)

    '''test QuantitativeFilterMixin here as well'''
    def test_validation(self):  
        '''test that this test returned correct compounds'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["descriptor"], OrdMolDescriptor.objects.get(heading="fs"), "Got a different result, {}, for the descriptor than expected".format(cleaned["descriptor"]))
        self.assertEqual(cleaned["operator"], "lt", "Got a different result, {}, for the operator than expected".format(cleaned["operator"]))
        self.assertEqual(int(cleaned["value"]), 3, "Got a different result, {}, for the value than expected".format(cleaned["value"]))

    def test_is_empty(self):
        self.form.is_valid() 
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got {} results intead of 1".format(fetched.count()))
        self.assertEqual(fetched[0].compound.abbrev, "water", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].compound.abbrev))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
class OrdinalFilterFormLE(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = OrdMolDescriptor.objects.get(heading="fs").pk
        self.formData["operator"] = "le" 
        self.formData["value"] = 3 

    def setUp(self): 
        self.setUpFormData()
        self.form = OrdinalFilterForm(data=self.formData)

    '''test QuantitativeFilterMixin here as well'''
    def test_validation(self):  
        '''test that this test returned correct compounds'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["descriptor"], OrdMolDescriptor.objects.get(heading="fs"), "Got a different result, {}, for the descriptor than expected".format(cleaned["descriptor"]))
        self.assertEqual(cleaned["operator"], "le", "Got a different result, {}, for the operator than expected".format(cleaned["operator"]))
        self.assertEqual(int(cleaned["value"]), 3, "Got a different result, {}, for the value than expected".format(cleaned["value"]))

    def test_is_empty(self):
        self.form.is_valid() 
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 3, "Got {} results intead of 3".format(fetched.count()))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
class CategoryFilterFormFun(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = CatMolDescriptor.objects.get(heading="arb").pk
        self.formData["value"] = str(CategoricalDescriptorPermittedValue.objects.get(value="fun").pk)

    def setUp(self): 
        self.setUpFormData()
        self.form = CategoryFilterForm(data=self.formData)        

    def test_validation(self):  
        '''test that this test returned correct compounds'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["descriptor"], CatMolDescriptor.objects.get(heading="arb"), "Got a different result, {}, for the descriptor than expected".format(cleaned["descriptor"]))
        self.assertEqual(int(cleaned["value"]), CategoricalDescriptorPermittedValue.objects.get(value="fun").pk, "Got {} instead of the expected value, {}".format(cleaned["value"], CategoricalDescriptorPermittedValue.objects.get(value="fun").pk))

    def test_is_empty(self):
        self.form.is_valid() 
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 4, "Got {} compounds instead of 4 compounds".format(fetched.count()))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
class CategoryFilterFormDull(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = CatMolDescriptor.objects.get(heading="arb").pk
        self.formData["value"] = CategoricalDescriptorPermittedValue.objects.get(value="dull").pk

    def setUp(self): 
        self.setUpFormData()
        self.form = CategoryFilterForm(data=self.formData)        

    def test_validation(self):  
        '''test that this test returned correct compounds'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["descriptor"], CatMolDescriptor.objects.get(heading="arb"), "Got a different result, {}, for the descriptor than expected".format(cleaned["descriptor"]))
        self.assertEqual(int(cleaned["value"]), CategoricalDescriptorPermittedValue.objects.get(value="dull").pk, "Got {}, rather than the expected result {}".format(cleaned["value"], CategoricalDescriptorPermittedValue.objects.get(value="dull").pk))

    def test_is_empty(self):
        self.form.is_valid() 
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 4, "Got {} instead of 4 compounds".format(fetched.count()))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
class BooleanFilterFormTest(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = BoolMolDescriptor.objects.get(heading="N?").pk
        self.formData["value"] = False

    def setUp(self): 
        self.setUpFormData()
        self.form = BooleanFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct compounds'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["descriptor"], BoolMolDescriptor.objects.get(heading="N?"), "Got a different result, {}, for the descriptor than expected".format(cleaned["descriptor"]))
        self.assertEqual(cleaned["value"], False, "Got a different result, {}, for the value than expected".format(cleaned["value"]))

    def test_is_empty(self):
        self.form.is_valid()
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got {} results instead of 1".format(fetched.count()))
        self.assertEqual(fetched[0].compound.abbrev, "water", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].compound.abbrev))


suite = unittest.TestSuite([
          loadTests(NumericFilterFormGT),
          loadTests(NumericFilterFormGE),
          loadTests(OrdinalFilterFormLT),
          loadTests(OrdinalFilterFormLE), 
          loadTests(CategoryFilterFormFun), 
          loadTests(CategoryFilterFormDull), 
          loadTests(BooleanFilterFormTest), 
          ])

if __name__=='__main__':
  runTests(suite)

