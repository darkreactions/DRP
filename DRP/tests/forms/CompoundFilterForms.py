#!/usr/bin/env python
'''Test classes for the compound administration form'''

import unittest 
from DRP.forms import filterForms  
from BaseFormTest import BaseFormTest
from DRP.models import LabGroup, ChemicalClass
from django.conf import settings
from django.contrib.auth.models import User
from DRP.tests.decorators import createsCompound, createsUser, joinsLabGroup, createsChemicalClass 

#@createsChemicalClass('Org', 'Organic Reagent')
#@createsCompound('EtOH', 682, 'Org', 'Narnia')

#create a filter request and see if the desired results are returned
#create a filter that is wrong and make sure it breaks (search for a compound that does not exist?)

loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
@createsChemicalClass('Org', 'Organic Reagent')
class CompoundFilterFormSucceed(BaseFormTest):
  '''This will test the fetch method of CompoundFilterForm class when it should succeed'''
  def setUpFormData(self):
   #custom, for Null, True, or False Boolean Field will be left out 

   self.formData["abbrevs"] = 
   self.formData["name"] =  
   self.formData["chemicalClasses"] = 
  def setUp(self):
  '''Creates a user, then a form'''
    self.setUpFormData()
    self.user = User.objects.get(username='Gamora') 

  def test_validation(self):
  '''create a filter request that should succeed, make sure it returns correct compounds'''  
    self.validationSucceeds() 



@createsUser('Gamora', 'pineapple_song')
@joinsLabGroup('Gamora', 'GalaxyGuardians')
class CompoundFilterFormFail(BaseFormTest)
  '''This will test the fetch method of CompoundFilterForm class when it should fail'''

  def setUp(self):
  '''Creates a user, then a form'''
    self.setUpFormData()
    self.user = User.objects.get(username='Gamora') 

  def test_validation(self):
  '''create a filter request that should fail'''  
    self.validationFails() 



class AdvancedCompoundFilterFormSuceed(BaseFormTest):

class AdvancedCompoundFilterFormFail(BaseFormTest): 



class NumericFilterFormSucceed(BaseFormTest):
  '''test QuantitativeFilterMixin here as well'''

class NumericFilterFormFail(BaseFormTest):
  '''test QuantitativeFilterMixin here as well'''



class OrdinalFilterFormSucceed(BaseFormTest):
  '''test QuantitativeFilterMixin here as well'''

class OrdinalFilterFormFail(BaseFormTest):
  '''test QuantitativeFilterMixin here as well'''



class CategoryFilterFormSucceed(BaseFormTest):
class CategoryFilterFormFail(BaseFormTest):



class BooleanFilterFormSucceed(BaseFormTest):
class BooleanFilterFormFail(BaseFormTest):


class CompoundFilterFormSetSucceed(BaseFormTest):
class CompoundFilterFormSetFail(BaseFormTest):


class AdvancedCompoundFilterFormSetSucceed(BaseFormTest):
class AdvancedCompoundFilterFormSetFail(BaseFormTest):
