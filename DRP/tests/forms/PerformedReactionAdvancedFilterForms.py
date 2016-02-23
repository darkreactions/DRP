#!/usr/bin/env python
'''Test classes for the performed reaction filter form'''

import unittest 
from DRP.forms import FilterForm, PerformedReactionFilterForm
from BaseFormTest import BaseFormTest
from DRP.models import LabGroup 
from django.conf import settings
from django.contrib.auth.models import User
from DRP.tests.decorators import createsUser, joinsLabGroup, loadsCompoundsFromCsv, _createsPerformedReaction
from DRP.tests.DRPTestCase import DRPTestCase, runTests
from DRP.models import CategoricalDescriptorPermittedValue
loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["128.2153"], {"Numeric": "128.2153", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class NumericFilterFormGT(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = NumMolDescriptor.objects.get(heading="mw").pk
        self.formData["operator"] = 'gt' 
        self.formData["value"] = 128.2153

    def setUp(self): 
        '''Creates a user, then a form'''
        self.setUpFormData()
        self.form = NumericFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid() 
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["descriptor"], NumMolDescriptor.objects.get(heading="mw"), "Got a different result, {}, for the descriptor than expected".format(cleaned["descriptor"]))
        self.assertEqual(cleaned["operator"], "gt", "Got a different result, {}, for the operator than expected".format(cleaned["descriptor"]))
        self.assertEqual(float(cleaned["value"]), 128.2153, "Got a different result, {}, for the value than expected".format(cleaned["value"]))

    def test_is_empty(self):
        self.form.is_valid() 
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got {} objects instead of 1".format(fetched.count()))
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["128.2153"], {"Numeric": "128.2153", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['hmta'], ["org"], ["140.2"], {"Numeric": "140.2", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "brt", "okr", False, True, True)
class NumericFilterFormGE(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = NumMolDescriptor.objects.get(heading="mw").pk
        self.formData["operator"] = 'ge' 
        self.formData["value"] = 128.2153

    def setUp(self): 
        '''Creates a user, then a form'''
        self.setUpFormData()
        self.form = NumericFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["descriptor"], NumMolDescriptor.objects.get(heading="mw"), "Got a different result, {}, for the descriptor than expected".format(cleaned["descriptor"]))
        self.assertEqual(cleaned["operator"], "ge", "Got a different result, {}, for the operator than expected".format(cleaned["operator"]))
        self.assertEqual(float(cleaned["value"]), 128.2153, "Got a different result,{}, for the value than expected".format(cleaned["value"]))

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
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "2", "Category": "fun"}, "abc", "2cd", False, True, False) 
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
        '''test that this test returned correct performed reactions'''
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
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
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
        '''test that this test returned correct performed reactions'''
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
        self.assertEqual(fetched.count(), 1, "Got {} results intead of 1".format(fetched.count()))
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class CategoryFilterFormFun(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = CatMolDescriptor.objects.get(heading="arb").pk
        self.formData["value"] = str(CategoricalDescriptorPermittedValue.objects.get(value="fun").pk)

    def setUp(self): 
        self.setUpFormData()
        self.form = CategoryFilterForm(data=self.formData)        

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
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
        self.assertEqual(fetched.count(), 1, "Got {} performed reactions instead of . performed reactions".format(fetched.count()))
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "dull"}, "abc", "2cd", False, True, False) 
class CategoryFilterFormDull(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = CatMolDescriptor.objects.get(heading="arb").pk
        self.formData["value"] = CategoricalDescriptorPermittedValue.objects.get(value="dull").pk

    def setUp(self): 
        self.setUpFormData()
        self.form = CategoryFilterForm(data=self.formData)        

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
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
        self.assertEqual(fetched.count(), 1, "Got {} instead of . performed reactions".format(fetched.count()))
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": False, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class BooleanFilterFormTest(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["descriptor"] = BoolMolDescriptor.objects.get(heading="N?").pk
        self.formData["value"] = False

    def setUp(self): 
        self.setUpFormData()
        self.form = BooleanFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
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
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class CompoundQuantityFilterFormTestCompound(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["compound"] =  "2-amep"

    def setUp(self): 
        self.setUpFormData()
        self.form = CompoundQuantityFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data

    def test_is_empty(self):
        self.form.is_valid()
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got {} results instead of 1".format(fetched.count()))
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class CompoundQuantityFilterFormTestRole(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["role"] =  "org"

    def setUp(self): 
        self.setUpFormData()
        self.form = CompoundQuantityFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data

    def test_is_empty(self):
        self.form.is_valid()
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got {} results instead of 1".format(fetched.count()))
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class CompoundQuantityFilterFormTestAmount(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["amount"] =  126.6

    def setUp(self): 
        self.setUpFormData()
        self.form = CompoundQuantityFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data

    def test_is_empty(self):
        self.form.is_valid()
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got {} results instead of 1".format(fetched.count()))
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False, "12/23/2013") 
class PerformedDateTimeTest(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["operator"] = "le" #on or before
        self.formData["date"] = "12/27/2013"

    def setUp(self): 
        self.setUpFormData()
        self.form = PerformedDateFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["operator"], "eq", "Got a different result, {}, for the descriptor than expected".format(cleaned["operator"]))

    def test_is_empty(self):
        self.form.is_valid()
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got {} results instead of 1".format(fetched.count()))
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test16.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False, "12/23/2013", "1/25/2014") 
class InsertedDateTimeTest(BaseFormTest):
    def setUpFormData(self):
        self.formData = {}
        self.formData["operator"] = "eq"
        self.formData["date"] = "1/25/2014"

    def setUp(self): 
        self.setUpFormData()
        self.form = InsertedDateTimeFilterForm(data=self.formData)

    def test_validation(self):  
        '''test that this test returned correct performed reactions'''
        self.validationSucceeds()     

    def test_clean(self):
        self.form.is_valid()
        cleaned = self.form.cleaned_data
        self.assertEqual(cleaned["operator"], "eq", "Got a different result, {}, for the descriptor than expected".format(cleaned["operator"]))
        self.assertEqual(cleaned["date"], "1/25/2014", "Got a different result, {}, for the descriptor than expected".format(cleaned["date"]))

    def test_is_empty(self):
        self.form.is_valid()
        self.assertFalse(self.form.is_empty())        

    def test_fetch(self): 
        self.form.is_valid()
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got {} results instead of 1".format(fetched.count()))
        self.assertEqual(fetched[0].performedreaction.reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].performedreaction.reference))



suite = unittest.TestSuite([
          loadTests(NumericFilterFormGT),
          loadTests(NumericFilterFormGE),
          loadTests(OrdinalFilterFormLT),
          loadTests(OrdinalFilterFormLE), 
          loadTests(CategoryFilterFormFun), 
          loadTests(CategoryFilterFormDull), 
          loadTests(BooleanFilterFormTest), 
          loadTests(CompoundQuantityFilterFormTestCompound),
          loadTests(CompoundQuantityFilterFormTestRole),
          loadTests(CompoundQuantityFilterFormTestAmount),
          loadTests(PerformedDateTimeTest),
          loadTests(InsertedDateTimeTest)
          ])

if __name__=='__main__':
  runTests(suite)

