#!/usr/bin/env python
'''Test classes for the compound filter form'''

import unittest 
from DRP.forms import FilterForm, PerformedReactionFilterForm
from BaseFormTest import BaseFormTest
from DRP.models import LabGroup 
from django.conf import settings
from django.contrib.auth.models import User
from DRP.tests.decorators import createsCompound, createsUser, joinsLabGroup, loadsCompoundsFromCsv, createsPerformedReactionVersion2
from DRP.tests.DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase
@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False)
class PerformedReactionFilterFormSucceed(BaseFormTest):
    '''This will test the PerformedReactionFilterForm class when it should succeed;
    Specifically, this test is creating a lab group with one performed reactin'''
    def setUpFormData(self):
        #custom, for Null, True, or False Boolean Field will be left out 
        self.formData = {}
        self.formData["reference"] = "abc" 
        self.formData["duplicateOf"] = "2cd"
        self.formData["legacyRecommendationFlag"] = False
        self.formData["valid"] = True
        self.formData["public"] = False 
        self.formData["labGroup"] = self.labgroup.pk 
        self.formData["js_active"] = False
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora') 
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationSucceeds() 

    def test_is_empty(self):
        self.assertTrue(self.form.is_valid(), str(self.form.errors))
        self.assertFalse(self.form.is_empty())

    def test_fetch(self):
        self.assertTrue(self.form.is_valid())
        fetched = self.form.fetch()
        self.assertEqual(fetched.count(), 1, "Got more than one result from the filter form")
        self.assertEqual(fetched[0].reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))
    
@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class PerformedReactionFilterFormDifferentLabGroupFail(BaseFormTest):
    '''This will test the PerformedReactionFilterForm when the labgroup searched for does not match the reaction searched for (i.e. searching for another labgroup's reactions), and the form should fail'''

    def setUpFormData(self):
        self.formData = {}
        self.formData["reference"] = "abc" 
        self.formData["duplicateOf"] = "2cd"
        self.formData["legacyRecommendationFlag"] = False
        self.formData["valid"] = True
        self.formData["public"] = False 
        self.formData["labGroup"]= self.labgroup.pk
        self.formData["js_active"]=False
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Starlord') 
        self.labgroup = LabGroup.objects.get(title='Terran')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''create a filter request that should fail, make sure it returns correct compounds'''  
        self.validationFails() 

    def test_is_empty(self):
        self.assertFalse(self.form.is_valid())
        self.assertFalse(self.form.is_empty())
        
@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'Alien')
@createsCompound("2-amep", "104820", "Amine", "Alien")
@createsUser('Starlord', 'gamora')
@joinsLabGroup('Starlord', 'Terran')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class PerformedReactionFilterFormWrongLabGroupFail(BaseFormTest):
    '''This will test the PerformedReactionFilterForm when the user tries to input a different labgroup in the form than the labgroup to which the user actually belongs'''

    def setUpFormData(self):
        self.formData = {}
        self.formData["reference"] = "abc" 
        self.formData["duplicateOf"] = "2cd"
        self.formData["legacyRecommendationFlag"] = False
        self.formData["valid"] = True
        self.formData["public"] = False 
        self.formData["labGroup"]= self.labgroup.pk
        self.formData["js_active"]=False
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Starlord') 
        self.labgroup = LabGroup.objects.get(title='Terran')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''Create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationFails() 

    def test_is_empty(self):
        self.assertFalse(self.form.is_valid())
        self.assertFalse(self.form.is_empty())
       
@createsUser('Starlord', 'gamora')
@joinsLabGroup('Starlord', 'Terran')
@loadsCompoundsFromCsv('Terran', 'compound_spread_test9.csv')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class PerformedReactionFilterFormEmptyFail(BaseFormTest):
    '''Test to make sure that PerformedReactionFilterForm fails with an empty form that does not specify the user/labgroup'''
    def setUpFormData(self): 
        self.formData = {} 
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Starlord')
        self.labgroup = LabGroup.objects.get(title='Terran')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''Create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationFails()  

    def test_is_empty(self):
        self.assertFalse(self.form.is_valid())
        self.assertTrue(self.form.is_empty())

####The following test one field at a time, (with the labgroup always filled out) to ensure each individual filter returns what is expected####

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class PerformedReactionFilterFormLabGroup(BaseFormTest):
    '''Test to make sure that PerformedReactionFilterForm succeeds/returns expected results with only the lab group specified'''
    def setUpFormData(self): 
        self.formData = {} 
        self.formData["labGroup"] = self.labgroup.pk
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora')
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationSucceeds() 

    def test_is_empty(self):
        self.assertTrue(self.form.is_valid())
        self.assertFalse(self.form.is_empty())

    def test_fetch(self):
        self.assertTrue(self.form.is_valid())
        fetched = self.form.fetch() 
        self.assertEqual(fetched.count(),1, "Expected 1 objects to be returned by form, number of objects returned was {}".format(fetched.count()))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class PerformedReactionFilterFormReference(BaseFormTest):
    '''Test that Performed Reaction Filter Form succeeds/returns expected results with only the reference name specified'''
    def setUpFormData(self): 
        self.formData = {} 
        self.formData["reference"] = "abc"
        self.formData["labGroup"] = self.labgroup.pk
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora') 
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''Create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationSucceeds() 

    def test_is_empty(self):
        self.assertTrue(self.form.is_valid())
        self.assertFalse(self.form.is_empty())

    def test_fetch(self):
        self.assertTrue(self.form.is_valid())
        fetched = self.form.fetch() 
        self.assertEqual(fetched.count(),1, "Expected 1 object to be returned by form, number of objects returned was {}".format(fetched.count()))
        self.assertEqual(fetched[0].reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class PerformedReactionFilterFormValid(BaseFormTest): 
    '''Test that the Performed Reaction Filter Form succeeds/returns expected results with only the 'valid' field specified'''
    def setUpFormData(self):
        self.formData = {}
        self.formData["valid"] = True 
        self.formData["labGroup"] = self.labgroup.pk
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora') 
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''Create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationSucceeds() 

    def test_is_empty(self):
        self.assertTrue(self.form.is_valid())
        self.assertFalse(self.form.is_empty())

    def test_fetch(self):
        self.assertTrue(self.form.is_valid())
        fetched = self.form.fetch() 
        self.assertEqual(fetched.count(),1, "Expected 1 object to be returned by form, number of objects returned was {}".format(fetched.count()))
        self.assertEqual(fetched[0].reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class PerformedReactionFilterFormDuplicateOf(BaseFormTest):
    '''Test that the Performed Reaction Filter Form succeeds/returns expected results with only the duplicateOf field specified'''

    def setUpFormData(self): 
        self.formData = {}
        self.formData["duplicateOf"]= "2cd"
        self.formData["labGroup"] = self.labgroup.pk
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora') 
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''Create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationSucceeds() 

    def test_is_empty(self):
        self.assertTrue(self.form.is_valid())
        self.assertFalse(self.form.is_empty())
    
    def test_fetch(self):
        self.assertTrue(self.form.is_valid())
        fetched = self.form.fetch() 
        self.assertEqual(fetched.count(),4, "Expected 4 objects to be returned by form, number of objects returned was {}".format(fetched.count()))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class PerformedReactionFilterFormLegacy(BaseFormTest): 
    '''Test that PerformedReaction Filter Form succeeds/returned expected results with only the legacyRecommendationFlag specified'''
    def setUpFormData(self):
        self.formData = {}
        self.formData["legacyRecommendationFlag"] = False
        self.formData["labGroup"] = self.labgroup.pk
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora') 
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''Create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationSucceeds() 

    def test_is_empty(self):
        self.assertTrue(self.form.is_valid())
        self.assertFalse(self.form.is_empty())

    def test_fetch(self):
        self.assertTrue(self.form.is_valid())
        fetched = self.form.fetch() 
        self.assertEqual(fetched.count(),1, "Expected 1 object to be returned by form, number of objects returned was {}".format(fetched.count()))
        self.assertEqual(fetched[0].reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@createsPerformedReactionVersion2('GalaxyGuardians', 'Gamora', ['2-amep'], ["Org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
class PerformedReactionFilterFormPublic(BaseFormTest):
    '''Test that Performed Reaction Filter Form succeeds/returned expected results with only the INCHI specified'''
    def setUpFormData(self):
        self.formData = {}
        self.formData["public"] = False
        self.formData["labGroup"] = self.labgroup.pk
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora') 
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.setUpFormData()
        self.form = PerformedReactionFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''Create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationSucceeds() 

    def test_is_empty(self):
        self.assertTrue(self.form.is_valid())
        self.assertFalse(self.form.is_empty())

    def test_fetch(self):
        self.assertTrue(self.form.is_valid())
        fetched = self.form.fetch() 
        self.assertEqual(fetched.count(),1, "Expected 1 object to be returned by form, number of objects returned was {}".format(fetched.count()))
        self.assertEqual(fetched[0].reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))



suite = unittest.TestSuite([
          loadTests(PerformedReactionFilterFormSucceed),
          loadTests(PerformedReactionFilterFormDifferentLabGroupFail),
          loadTests(PerformedReactionFilterFormWrongLabGroupFail),
          loadTests(PerformedReactionFilterFormEmptyFail), 
          loadTests(PerformedReactionFilterFormLabGroup),
          loadTests(PerformedReactionFilterFormReference),
          loadTests(PerformedReactionFilterFormValid),
          loadTests(PerformedReactionFilterFormDuplicateOf),
          loadTests(PerformedReactionFilterFormLegacy),
          loadTests(PerformedReactionFilterFormPublic)
          ])

if __name__=='__main__':
  runTests(suite)

