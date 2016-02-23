#!/usr/bin/env python
'''Test classes for the compound filter form'''

import unittest 
from DRP.forms import FilterForm, PerformedReactionFilterFormSet
from BaseFormTest import BaseFormTest
from DRP.models import LabGroup
from django.conf import settings
from django.contrib.auth.models import User
from DRP.tests.decorators import createsCompound, createsUser, joinsLabGroup, loadsCompoundsFromCsv, _createsPerformedReaction
from DRP.tests.DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase
@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['2-amep'], ["org"], ["141.5"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "abc", "2cd", False, True, False) 
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['hmta'], ["org"], ["126.6"], {"Numeric": "77", "Boolean": True, "Ordinal": "3", "Category": "fun"}, "brt", "okr", False, True, True)
@_createsPerformedReaction('GalaxyGuardians', 'Gamora', ['sov'], ["org"], ["75"], {"Numeric": "77", "Boolean": False, "Ordinal": "3", "Category": "fun"}, "sup", "hey", True, True, True)
class PerformedReactionFilterSetSucceed(BaseFormTest):
    '''Three different forms searching for different compounds, should all succeed'''
    def setUpFormData(self):
        self.formData = {}
        self.formData["form-TOTAL_FORMS"] = 3
        self.formData["form-INITIAL_FORMS"] = 0
        self.formData["form-MAX_NUM_FORMS"] = 10000
        #first form
        self.formData["form-0-reference"] = "abc"
        self.formData["form-0-duplicateOf"]= "2cd"
        self.formData["form-0-legacyRecommendationFlag"] = False
        self.formData["form-0-valid"] =	True
        self.formData["form-0-public"] = False
        self.formData["form-0-labGroup"] = self.labgroup.pk 
        self.formData["form-0-js_active"] = False
	    #second form 
        self.formData["form-0-reference"] = "brt"
        self.formData["form-0-duplicateOf"]= "okr"
        self.formData["form-0-legacyRecommendationFlag"] = False
        self.formData["form-0-valid"] =	True
        self.formData["form-0-public"] = True
        self.formData["form-1-labGroup"] = self.labgroup.pk 
        self.formData["form-1-js_active"] = False
        #third form 	
        self.formData["form-0-reference"] = "sup"
        self.formData["form-0-duplicateOf"]= "hey"
        self.formData["form-0-legacyRecommendationFlag"] = True
        self.formData["form-0-valid"] =	True
        self.formData["form-0-public"] = True
        self.formData["form-2-labGroup"] = self.labgroup.pk 
        self.formData["form-2-js_active"] = False
        
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora')
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.setUpFormData()
        self.form = PerformedReactionFilterFormSet(user=self.user, labGroup=self.labgroup, data=self.formData)

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
        self.assertEqual(fetched[0].reference, "abc", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))
        self.assertEqual(fetched[1].reference, "brt", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))
        self.assertEqual(fetched[2].reference, "sup", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))
        self.assertEqual(fetched.count(), 3, "{} results found.".format(fetched.count()))
        

suite = unittest.TestSuite([
          loadTests(PerformedReactionFilterSetSucceed), 
          ])

if __name__=='__main__':
  runTests(suite)

