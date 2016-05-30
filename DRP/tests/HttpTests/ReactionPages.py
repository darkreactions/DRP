#!/usr/bin/env python

from django.conf import settings
from HttpTest import GetHttpTest, PostHttpTest, GetHttpSessionTest, PostHttpSessionTest
from HttpTest import  redirectionMixinFactory, logsInAs, usesCsrf
from HttpTest import  choosesLabGroup 
from DRP.tests.decorators import joinsLabGroup, createsChemicalClass, signsExampleLicense
from DRP.tests.decorators import createsUser, createsCompound, createsPerformedReaction
from DRP.tests import runTests
from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
import requests
import unittest
from DRP.models import LabGroup, CompoundRole, Compound, PerformedReaction

newReactionUrl = reverse('newReaction')
reactionBaseCodes = ['ea5108a2-0b88-482d-90c5-ea492fd8134e', '823d22e7-3337-4292-aa67-13f748b2aa65', '1f47e7ab-1900-4683-ba1c-63330ec2f71a']

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
class GetReactionCreate(GetHttpSessionTest):

    url=newReactionUrl
    testCodes=['6813e404-f1a3-48ed-ae97-600a41cf63cb','2758c44c-b7e2-440a-a617-36d9d730bc93'] + reactionBaseCodes


@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@usesCsrf
class PostReactionCreateValid(PostHttpSessionTest, redirectionMixinFactory(1)): # effectively tests GET for add_reactants view

    url=newReactionUrl
    testCodes=['008d2580-5be2-4112-8297-a9e53490bb6d', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06'] + reactionBaseCodes
    _payload={'reference':'turkish_delight'}

    def setUp(self):
        self.payload['labGroup']=LabGroup.objects.get(title='Narnia').id
        super(PostReactionCreate, self).setUp()

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@usesCsrf
class PostReactionCreateInvalid(PostHttpSessionTest):

    url=newReactionUrl
    testCodes=['6813e404-f1a3-48ed-ae97-600a41cf63cb','2758c44c-b7e2-440a-a617-36d9d730bc93'] + reactionBaseCodes
    _payload={'reference':''}

    def setUp(self):
        self.payload['labGroup']=LabGroup.objects.get(title='Narnia').id
        super(PostReactionCreate, self).setUp()

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@createsPerfReaction('narnia', 'Aslan', 'turkish_delight')
class GetReactionEdit(GetHttpSessionTest):

    testCodes=['7b3b6668-981a-4a11-8dc4-23107187de93', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06', '2758c44c-b7e2-440a-a617-36d9d730bc93'] + reactionBaseCodes

    def __init__(*args, **kwargs):
        self.url = reverse('editReaction', kwargs={'rxn_id':PreformedReaction.objects.get(reference='turkish_delight', labGroup__title='narnia')})
        super(GetReactionEdit, self).__init__(*args, **kwargs)


@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@createsPerfReaction('narnia', 'Aslan', 'turkish_delight')
@createsOrdRxnDescriptor('deliciousness', 0, 4)
class GetReactionEdit2(GetHttpSessionTest):

    url=editReactionUrl
    testCodes=['7b3b6668-981a-4a11-8dc4-23107187de93', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06', '634d88bb-9289-448b-a3dc-548ff4c6cda1', '2758c44c-b7e2-440a-a617-36d9d730bc93'] + reactionBaseCodes

@logsInAs('Aslan', 'old_magic')
@createsUser('WhiteQueen', 'New Magic')
@signsExampleLicense('Aslan')
@signsExampleLicense('WhiteQueen')
@joinsLabGroup('Aslan', 'narnia')
@joinsLabGroup('WhiteQueen', 'WhiteQueensArmy')
@createsPerfReaction('narnia', 'Aslan', 'turkish_delight')
@createsPerfReaction('WhiteQueensArmy', 'WhiteQueen', 'turkish_delight')
class GetSomeoneElsesReactionEdit(GetHttpSessionTest):
    '''Attempts to grab someone else\'s reaction, which should fail with a 404.''''

    status=404

    def __init__(*args, **kwargs):
        self.url = reverse('editReaction', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight', labGroup__title='WhiteQueensArmy').id})
        super(GetReactionEdit, self).__init__(*args, **kwargs)

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
class GetNonexistentReactionEdit(GetHttpSessionTest):

    status=404

    def __init__(*args, **kwargs):
        self.url = reverse('editReaction', kwargs={'rxn_id':4})
        super(GetReactionEdit, self).__init__(*args, **kwargs)

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@createsPerfReaction('narnia', 'Aslan', 'turkish_delight')
@usesCsrf
class PostReactionEditValid(PostHttpSessionTest):

    testCodes=['7b3b6668-981a-4a11-8dc4-23107187de93', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06', '2758c44c-b7e2-440a-a617-36d9d730bc93', 'd96fc7a1-69cf-44ac-975d-a67f9e2c74d0'] + reactionBaseCodes

    def setUp(self):
        self.url = reverse('editReaction', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight', labGroup__title='narnia').id})
        self.reaction = PerformedReaction.objects.get('turkish_delight', labGroup__title='narnia')
        self.payload.update({'notes':'this reaction has been edited', 'reference':self.reaction.reference, 'labGroup':self.reaction.labGroup.id, 'performedBy':self.reaction.user.id, 'performedDateTime':self.reaction.performedDateTime, 'valid':self.reaction.valid, 'duplicateOf':""})
    super(PostReactionEditValid, self).setUp()

    def test_edit(self):
        self.assertEqual(self.reaction.notes,'this reaction has been edited')

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@createsPerfReaction('narnia', 'Aslan', 'turkish_delight')
@usesCsrf
class PostReactionEditInvalid(PostHttpSessionTest):

    testCodes=['7b3b6668-981a-4a11-8dc4-23107187de93', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06', '2758c44c-b7e2-440a-a617-36d9d730bc93'] + reactionBaseCodes

    def setUp(self):
        self.url = reverse('editReaction', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight', labGroup__title='narnia').id})
        self.reaction = PerformedReaction.objects.get('turkish_delight', labGroup__title='narnia')
        self.payload.update({'notes':'this reaction has been edited', 'labGroup':self.reaction.labGroup.id, 'performedBy':self.reaction.user.id, 'performedDateTime':self.reaction.performedDateTime, 'valid':self.reaction.valid, 'duplicateOf':""})
    super(PostReactionEditValid, self).setUp()
        
    def test_edit(self):
        self.assertNotEqual(self.reaction.notes,'this reaction has been edited')

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@createsPerfReaction('narnia', 'Aslan', 'turkish_delight')
@usesCsrf
class DeletePerformedReaction(PostHttpSessionTest, redirectionMixinFactory(1)):

    testCodes = ['ba960469-5fad-4142-a0a1-a10b37e9432e'] + reactionBaseCodes
    url = reverse('deleteReaction')

    def setUp(self):
        self.reaction = PerformedReaction.objects.get('turkish_delight', labGroup__title='narnia')
        self.payload['id'] = self.reaction.id
        super(DeletePerformedReaction, self).setUp()

    def tests_deleted(self):
        self.assertFalse(PerformedReaction.objects.filter('turkish_delight', labGroup__title='narnia').exists())

@logsInAs('Aslan', 'old_magic')
@createsUser('WhiteQueen', 'New Magic')
@signsExampleLicense('Aslan')
@signsExampleLicense('WhiteQueen')
@joinsLabGroup('Aslan', 'narnia')
@joinsLabGroup('WhiteQueen', 'WhiteQueensArmy')
@createsPerfReaction('narnia', 'Aslan', 'turkish_delight')
@createsPerfReaction('WhiteQueensArmy', 'WhiteQueen', 'turkish_delight')
@usesCsrf
class DeleteSomeoneElsesReaction(PostHttpSessionTest, redirectionMixinFactory(1))
    testCodes = ['ba960469-5fad-4142-a0a1-a10b37e9432e'] + reactionBaseCodes
    url = reverse('deleteReaction')

    def setUp(self):
        self.reaction = PerformedReaction.objects.get('turkish_delight', labGroup__title='WhiteQueensArmy')
        self.payload['id'] = self.reaction.id
        super(DeleteSomeoneElsesReaction, self).setUp()

    def tests_deleted(self):
        self.assertTrue(PerformedReaction.objects.filter('turkish_delight', labGroup__title='narnia').exists())

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@usesCsrf
class DeleteNonexistentReaction(PostHttpSessionTest, redirectionMixinFactory(1)):
    url = reverse('deleteReaction')

    testCodes = ['ba960469-5fad-4142-a0a1-a10b37e9432e'] + reactionBaseCodes

    def setUp(self):
        self.payload['id'] = 5
        super(DeleteNonExistentReaction, self).setUp()
    
@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@createsPerfReaction('narnia', 'Aslan', 'turkish_delight')
@usesCsrf
class InvalidatePerformedReaction(PostHttpSessionTest, redirectionMixinFactory(1)):

    url = reverse('invalidateReaction')
    testCodes = ['ba960469-5fad-4142-a0a1-a10b37e9432e'] + reactionBaseCodes

    def setUp(self):
        self.reaction = PerformedReaction.objects.get('turkish_delight', labGroup__title='narnia')
        self.payload['id'] = self.reaction.id
        super(InvalidatePerformedReaction, self).setUp()

    def tests_deleted(self):
        self.assertFalse(PerformedReaction.objects.get('turkish_delight', labGroup__title='narnia').valid)

@logsInAs('Aslan', 'old_magic')
@createsUser('WhiteQueen', 'New Magic')
@signsExampleLicense('Aslan')
@signsExampleLicense('WhiteQueen')
@joinsLabGroup('Aslan', 'narnia')
@joinsLabGroup('WhiteQueen', 'WhiteQueensArmy')
@createsPerfReaction('narnia', 'Aslan', 'turkish_delight')
@createsPerfReaction('WhiteQueensArmy', 'WhiteQueen', 'turkish_delight')
@usesCsrf
class InvalidateSomeoneElsesReaction(PostHttpSessionTest, redirectionMixinFactory(1))
    url = reverse('invalidateReaction')
    testCodes = ['ba960469-5fad-4142-a0a1-a10b37e9432e'] + reactionBaseCodes

    def setUp(self):
        self.reaction = PerformedReaction.objects.get('turkish_delight', labGroup__title='WhiteQueensArmy')
        self.payload['id'] = self.reaction.id
        super(InvalidateSomeoneElsesReaction, self).setUp()

    def tests_deleted(self):
        self.assertTrue(PerformedReaction.objects.get('turkish_delight', labGroup__title='narnia').valid)

@logsInAs('Aslan', 'old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@usesCsrf
class InvalidateNonExistentReaction(PostHttpSessionTest, redirectionMixinFactory(1)):

    url = reverse('invalidateReaction')
    testCodes = ['ba960469-5fad-4142-a0a1-a10b37e9432e'] + reactionBaseCodes

    def setUp(self):
        self.payload['id'] = 5
        super(InvalidateNonExistentReaction, self).setUp()

@logsInAs('Aslan, old_magic')
@signsExampleLicense('Aslan')
@joinsLabGorup('Aslan', 'narnia')
@loadsCompoundsFromCsv('narnia', 'compound_spread_test1.csv')
@createsPerformedReaction('Aslan', 'narnia', 'turkish_delight')
@createsCompoundRole('Org', 'Organic')
@createsCompoundRole('inOrg', 'inOrganic')
@usesCsrf
class PostReactantAddCreatingValid(PostHttpSessionTest, redirectionMixinFactory(4)): #we expect 4 redirections because we have initialised no manual reaction descriptors
    url=reverse('addCompoundDetails')
    testCodes=['7b3b6668-981a-4a11-8dc4-23107187de93'] + reactionBaseCodes
    _params={'creating':True}
    _payload={'amount':1}

    def setUp(self):
        self.payload['form-0-role']=CompoundRole.objects.get(label='Org').id
        self.payload['form-0-compound']=Compound.objects.get(abbrev='2-amep').id
        self.payload['form-0-amount']='22'
        super(PostReactatAddCreating, self).setUp()

@logsInAs('Aslan, old_magic')
@signsExampleLicense('Aslan')
@joinsLabGorup('Aslan', 'narnia')
@loadsCompoundsFromCsv('narnia', 'compound_spread_test1.csv')
@createsPerformedReaction('Aslan', 'narnia', 'turkish_delight')
@createsCompoundRole('Org', 'Organic')
@createsCompoundRole('inOrg', 'inOrganic')
@createsOrdRxnDescriptor('deliciousness', 0, 4)
@usesCsrf
class PostReactantAddCreatingValid2(PostHttpSessionTest, redirectionMixinFactory(2)): #we expect 4 redirections because we have initialised no manual reaction descriptors
    testCodes=['9fa2cfb6-aabe-40f7-80ea-4ecbcf8c0bda', '634d88bb-9289-448b-a3dc-548ff4c6cda1'] + reactionBaseCodes
    _params={'creating':True}
    _payload={'amount':1}

    def setUp(self):
        self.url=reverse('addCompoundDetails', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight')})
        self.payload['form-0-role']=CompoundRole.objects.get(label='Org').id
        self.payload['form-0-compound']=Compound.objects.get(abbrev='2-amep').id
        self.payload['form-0-amount']='22'
        super(PostReactantAddCreatingValid2, self).setUp()

@logsInAs('Aslan, old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@loadsCompoundsFromCsv('narnia', 'compound_spread_test1.csv')
@createsPerformedReaction('Aslan', 'narnia', 'turkish_delight')
@createsCompoundRole('Org', 'Organic')
@createsCompoundRole('inOrg', 'inOrganic')
@usesCsrf
class PostReactantAddCreatingInValid(PostHttpSessionTest):
    _params={'creating':True}
    testCodes=['008d2580-5be2-4112-8297-a9e53490bb6d', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06'] + reactionBaseCodes
    _payload={'amount':1}

    def setUp(self):
        self.url=reverse('addCompoundDetails', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight')})
        self.payload['form-0-role']=CompoundRole.objects.get(label='Org').id
        self.payload['form-0-compound']=Compound.objects.get(abbrev='2-amep').id
        self.payload['form-0-amount']='-1'
        super(PostReactantAddCreating, self).setUp()

@logsInAs('Aslan, old_magic')
@signsExampleLicense('Aslan')
@joinsLabGroup('Aslan', 'narnia')
@loadsCompoundsFromCsv('narnia', 'compound_spread_test1.csv')
@createsPerformedReaction('Aslan', 'narnia', 'turkish_delight')
@createsCompoundRole('Org', 'Organic')
@createsCompoundRole('inOrg', 'inOrganic')
@createsCompoundQuantity('turkish_delight', '2-amep', 'Org', '2.12')
@usesCsrf
class PostReactantEditingValid(PostHttpSessionTest, redirectionMixinFactory(1)):  #this is the same view so the invalid case is covered
    testCodes=['008d2580-5be2-4112-8297-a9e53490bb6d', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06'] + reactionBaseCodes
    _payload={'amount':1}

    def setUp(self):
        self.url=reverse('addCompoundDetails', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight')})
        self.quantity = CompoundQuantity.objects.filter(reaction__reference='turkish_delight', reaction__labGroup__title='narnia', compound_abbrev='2-amep')
        self.payload['form-0-id']=self.quantity.id
        self.payload['form-0-role']=CompoundRole.objects.get(label='Org').id
        self.payload['form-0-compound']=Compound.objects.get(abbrev='2-amep').id
        self.payload['form-0-amount']='15'
        super(PostReactantEditing, self).setUp()

    def test_edit(self):
        self.assertEqual(self.quantity.amount=15)


@logsInAs('Aslan, old_magic')
@signsExampleLicense('Aslan')
@joinsLabGorup('Aslan', 'narnia')
@loadsCompoundsFromCsv('narnia', 'compound_spread_test1.csv')
@createsPerformedReaction('Aslan', 'narnia', 'turkish_delight')
@createsCompoundRole('Org', 'Organic')
@createsCompoundRole('inOrg', 'inOrganic')
@createsCompoundQuantity('turkish_delight', '2-amep', 'Org', '2.12')
@usesCsrf
class PostReactantDeleteValid(PostHttpSessionTest, redirectionMixinFactory(1)): #we expect 4 redirections because we have initialised no manual reaction descriptors
    testCodes=['7b3b6668-981a-4a11-8dc4-23107187de93'] + reactionBaseCodes

    def setUp(self):
        self.url=reverse('addCompoundDetails', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight')})
        self.payload['form-0-role']=CompoundRole.objects.get(label='Org').id
        self.payload['form-0-compound']=Compound.objects.get(abbrev='2-amep').id
        self.payload['form-0-amount']='22'
        self.payload['form-0-delete']=''
        self.payload['form-0-id'] = CompoundQuantity.objects.get(compound_abbrev='2-amep', reaction__abbrev='turkish_delight')
        super(PostReactantDeleteValid, self).setUp()

    def test_delete(self):
        self.assertFalse(CompoundQuantity.objects.filter(compound_abbrev='2-amep', reaction__abbrev='turkish_delight').exists())
        
@logsInAs('Aslan, old_magic')
@signsExampleLicense('Aslan')
@joinsLabGorup('Aslan', 'narnia')
@loadsCompoundsFromCsv('narnia', 'compound_spread_test1.csv')
@createsPerformedReaction('Aslan', 'narnia', 'turkish_delight')
@createsOrdRxnDescriptor('deliciousness', 0, 4)
class CreateReactionDescValCreating(PostHttpSessionTest, redirectionMixinFactory(3)):
    testCodes=['008d2580-5be2-4112-8297-a9e53490bb6d', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06'] + reactionBaseCodes
    _params={'creating':True}
    
    def setUp(self):
        self.url=reverse('createOrdDescVals', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight')})
        self.payload['form-0-descriptor'] = OrdRxnDescriptor.objects.get(heading='deliciousness')
        self.payload['form-0-value']='0'
        super(CreateReactionDescValCreating, self).setUp()

    def test_created(self):
        self.assertTrue(OrdRxnDescriptorValue.objects.all().exists())
    
@logsInAs('Aslan, old_magic')
@signsExampleLicense('Aslan')
@joinsLabGorup('Aslan', 'narnia')
@loadsCompoundsFromCsv('narnia', 'compound_spread_test1.csv')
@createsPerformedReaction('Aslan', 'narnia', 'turkish_delight')
@createsOrdRxnDescriptor('deliciousness', 0, 4)
class CreateReactionDescValEditing(PostHttpSessionTest, redirectionMixinFactory(1)): #no edit test required- it uses the same functionality
    testCodes=['008d2580-5be2-4112-8297-a9e53490bb6d', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06'] + reactionBaseCodes
    
    def setUp(self):
        self.url=reverse('createOrdDescVals', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight')})
        self.payload['form-0-descriptor'] = OrdRxnDescriptor.objects.get(heading='deliciousness')
        self.payload['form-0-value']='0'
        super(CreateReactionDescValCreating, self).setUp()

    def test_created(self):
        self.assertTrue(OrdRxnDescriptorValue.objects.all().exists())

@logsInAs('Aslan, old_magic')
@signsExampleLicense('Aslan')
@joinsLabGorup('Aslan', 'narnia')
@loadsCompoundsFromCsv('narnia', 'compound_spread_test1.csv')
@createsPerformedReaction('Aslan', 'narnia', 'turkish_delight')
@createsOrdRxnDescriptor('deliciousness', 0, 4)
@createsOrdRxnDescriptorValue('narnia', 'turkish_delight', 'deliciousness', 3)
class DeleteReactionDescVal(PostHttpSessionTest, redirectionMixinFactory(1)):
    testCodes=['008d2580-5be2-4112-8297-a9e53490bb6d', 'dc1d5961-a9e7-44d8-8441-5b8402a01c06'] + reactionBaseCodes
    
    def setUp(self):
        self.url=reverse('createOrdDescVals', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight')})
        self.payload['form-0-descriptor'] = OrdRxnDescriptor.objects.get(heading='deliciousness')
        self.payload['form-0-value']='0'
        self.payload['form-0-DELETE']=''
        super(CreateReactionDescValCreating, self).setUp()

    def test_deleted(self):
        self.assertFalse(OrdRxnDescriptorValue.objects.all().exists())


class EditReactionDescValInvalid(PostHttpSessionTest):
@signsExampleLicense('Aslan')
@joinsLabGorup('Aslan', 'narnia')
@loadsCompoundsFromCsv('narnia', 'compound_spread_test1.csv')
@createsPerformedReaction('Aslan', 'narnia', 'turkish_delight')
@createsOrdRxnDescriptor('deliciousness', 0, 4)
@createsOrdRxnDescriptorValue('narnia', 'turkish_delight', 'deliciousness', 3)
class EditReactionDescValInvalid(PostHttpSessionTest):
    testCodes=['9fa2cfb6-aabe-40f7-80ea-4ecbcf8c0bda', '634d88bb-9289-448b-a3dc-548ff4c6cda1'] + reactionBaseCodes
    
    def setUp(self):
        self.url=reverse('createOrdDescVals', kwargs={'rxn_id':PerformedReaction.objects.get(reference='turkish_delight')})
        self.payload['form-0-descriptor'] = OrdRxnDescriptor.objects.get(heading='deliciousness')
        self.payload['form-0-value']='-1'
        super(CreateReactionDescValCreating, self).setUp()

    def test_created(self):
        self.assertTrue(OrdRxnDescriptorValue.objects.all().exists())

