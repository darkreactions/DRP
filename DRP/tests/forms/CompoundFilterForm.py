#!/usr/bin/env python
'''Test classes for the compound administration form'''

import unittest 
from DRP.forms import FilterForm, CompoundFilterForm 
from BaseFormTest import BaseFormTest
from DRP.models import LabGroup, ChemicalClass
from django.conf import settings
from django.contrib.auth.models import User
from DRP.tests.decorators import createsCompound, createsUser, joinsLabGroup, createsChemicalClass, loadsCompoundsFromCsv
from DRP.tests.DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase

#@createsChemicalClass('Org', 'Organic Reagent')
#@createsCompound('EtOH', 682, 'Org', 'Narnia')

#create a filter request and see if the desired results are returned
#create a filter that is wrong and make sure it breaks (search for a compound that does not exist?)

loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@createsChemicalClass("Amine", "descr")
@loadsCompoundsFromCsv('GalaxyGuardians', 'compound_spread_test1.csv')
class CompoundFilterFormSucceed(BaseFormTest):
    '''This will test the CompoundFilterForm class when it should succeed;
    Specifically, this test is creating a lab group with some chemical compounds, and 
    searching for one of those compounds (meaning it should definitely be returned)'''
    def setUpFormData(self):
        #custom, for Null, True, or False Boolean Field will be left out 
        self.formData = {}
        self.formData["abbrevs"] = "2-amep"
        self.formData["name"] =	"2-aminomethyl-1-ethylpyrrolidine"
        self.formData["chemicalClasses"]= [self.chemicalClass.pk]
        self.formData["CSID"] = "104820"
        self.formData["INCHI"] = "InChI=1S/C7H16N2/c1-2-9-5-3-4-7(9)6-8/h7H,2-6,8H2,1H3"
        self.formData["smiles"] = "CCN1CCCC1CN"
        self.formData["labGroup"] = self.labgroup.pk 
        self.formData["js_active"] = False
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Gamora') 
        self.labgroup = LabGroup.objects.get(title='GalaxyGuardians')
        self.chemicalClass = ChemicalClass.objects.get(label='Amine')
        self.setUpFormData()
        self.form = CompoundFilterForm(self.user, self.labgroup, self.formData) 

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
        self.assertEqual(fetched[0].abbrev, "2-amep", "The result returned did not have the expected abbrev attr. Attribute given was {}".format(fetched[0].abbrev))
    

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'Alien')
@createsChemicalClass("Amine", "descr")
@createsCompound("2-amep", "104820", "Amine", "Alien")
@createsUser('Starlord', 'gamora')
@joinsLabGroup('Starlord', 'Terran')
class CompoundFilterFormDifferentLabGroupFail(BaseFormTest):
    '''This will test the CompoundFilterForm when the labgroup searched for does not match the compound searched for (i.e. searching for another labgroup's compounds), and the form should fail'''
    def setUpFormData(self):
        self.formData={}
        self.formData["abbrevs"]="2-amep"
        self.formData["name"]="2-aminomethyl-1-ethylpyrrolidine"
        self.formData["chemicalClasses"]= [self.chemicalClass.pk]
        self.formData["CSID"]="104820"
        self.formData["INCHI"]="InChI=1S/C7H16N2/c1-2-9-5-3-4-7(9)6-8/h7H,2-6,8H2,1H3"
        self.formData["smiles"]="CCN1CCCC1CN"
        self.formData["labGroup"]= self.labgroup.pk
        self.formData["js_active"]=False
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Starlord') 
        self.labgroup = LabGroup.objects.get(title='Terran')
        self.chemicalClass = ChemicalClass.objects.get(label='Amine')
        self.setUpFormData()
        self.form = CompoundFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationFails() 

    def test_is_empty(self):
        self.assertFalse(self.form.is_valid())
        self.assertFalse(self.form.is_empty())
        
@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'Alien')
@createsChemicalClass("Amine", "descr") 
@createsCompound("2-amep", "104820", "Amine", "Alien")
@createsUser('Starlord', 'gamora')
@joinsLabGroup('Starlord', 'Terran')
class CompoundFilterFormWrongLabGroupFail(BaseFormTest):
    '''This will test the CompoundFilterForm when the user tries to input a different labgroup in the form than the labgroup to which the user actually belongs'''
    def setUpFormData(self): 
        self.formData = {} 
        self.formData["abbrevs"] = "2-amep"
        self.formData["name"] =  "2-aminomethyl-1-ethylpyrrolidine"
        self.formData["chemicalClasses"]= [self.chemicalClass.pk]
        self.formData["CSID"] = "104820"
        self.formData["INCHI"] = "InChI=1S/C7H16N2/c1-2-9-5-3-4-7(9)6-8/h7H,2-6,8H2,1H3"
        self.formData["smiles"] = "CCN1CCCC1CN"
        self.formData["labGroup"] = self.labgroup.pk
        self.formData["js_active"] = False
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Starlord') 
        self.labgroup = LabGroup.objects.get(title='Terran')
        self.chemicalClass = ChemicalClass.objects.get(label='Amine')
        self.setUpFormData()
        self.form = CompoundFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''Create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationFails() 

    def test_is_empty(self):
        self.assertFalse(self.form.is_valid())
        self.assertFalse(self.form.is_empty())

  
@createsUser('Starlord', 'gamora')
@joinsLabGroup('Starlord', 'Terran')
@loadsCompoundsFromCsv('Terran', 'compound_spread_test9.csv')
class CompoundFilterFormEmptySucceed(BaseFormTest):
    '''Test to make sure that CompoundFilterForm succeeds with an empty form'''
    def setUpFormData(self): 
        self.formData = {} 
        self.formData["labGroup"] = self.labgroup.pk
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Starlord')
        self.labgroup = LabGroup.objects.get(title='Terran')
        self.chemicalClass = ChemicalClass.objects.get(label='Amine')
        self.setUpFormData()
        self.form = CompoundFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationSucceeds() 

    def test_is_empty(self):
        self.assertTrue(self.form.is_valid())
        self.assertTrue(self.form.is_empty())

    def test_fetch(self):
        self.assertTrue(self.form.is_valid())
        fetched = self.form.fetch() 
        self.assertEqual(fetched.count(),8, "Expected 8 objects to be returned by form, number of objects returned was {}".format(fetched.count()))
        
@createsUser('Starlord', 'gamora')
@joinsLabGroup('Starlord', 'Terran')
@loadsCompoundsFromCsv('Terran', 'compound_spread_test9.csv')
class CompoundFilterFormEmptyFail(BaseFormTest):
    '''Test to make sure that CompoundFilterForm fails with an empty form that does not specify the user/labgroup'''
    def setUpFormData(self): 
        self.formData = {} 
	
    def setUp(self):
        '''Creates a user, then a form'''
        self.user = User.objects.get(username='Starlord')
        self.labgroup = LabGroup.objects.get(title='Terran')
        self.chemicalClass = ChemicalClass.objects.get(label='Amine')
        self.setUpFormData()
        self.form = CompoundFilterForm(self.user, self.labgroup, self.formData) 

    def test_validation(self):
        '''Create a filter request that should succeed, make sure it returns correct compounds'''  
        self.validationFails()  

    def test_is_empty(self):
        self.assertFalse(self.form.is_valid())
        self.assertTrue(self.form.is_empty())

suite = unittest.TestSuite([
          loadTests(CompoundFilterFormSucceed),
          loadTests(CompoundFilterFormDifferentLabGroupFail),
          loadTests(CompoundFilterFormWrongLabGroupFail),
          loadTests(CompoundFilterFormEmptySucceed),
          loadTests(CompoundFilterFormEmptyFail), 
          ])

if __name__=='__main__':
  runTests(suite)

