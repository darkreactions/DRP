from django.core.management.base import BaseCommand
from DRP.tests import suite, runTests
from django.conf import settings
from os import path
import csv
from django.contrib.auth.models import User
from DRP.models import LabGroup, Compound, PerformedReaction 
from DRP.models import ChemicalClass, NumRxnDescriptor
from DRP.models import BoolRxnDescriptor, OrdRxnDescriptor
from DRP.models import NumRxnDescriptorValue, BoolRxnDescriptorValue
from DRP.models import OrdRxnDescriptorValue, NumMolDescriptorValue 
from chemspipy import ChemSpider

outcomeDescriptor = OrdRxnDescriptor.objects.get_or_create(
    heading='crystallisation_outcome',
    name='Four class crystallisation outcome',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0',
    maximum=4,
    minimum=1
    )

outcomeBooleanDescriptor = BoolRxnDescriptor.objects.get_or_create(
    heading='boolean_crystallisation_outcome',
    name='Two class crystallisation outcome',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'    
    )

purityDescriptor = OrdRxnDescriptor.objects.get_or_create(
    heading='crystallisation_purity_outcome',
    name='Two class crystallisation purity',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0',
    maximum=2
    minimum=1
    )

temperatureDescriptor = NumRxnDescriptor.objects.get_or_create(
    heading='reaction_temperature',
    name='Reaction temperature in degrees C',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0' 
    )

timeDescriptor = NumRxnDescriptor.objects.get_or_create(
    heading='reaction_time',
    name='Reaction time in minutes',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'  
    )

pHDescriptor = NumRxnDescriptor.objects.get_or_create(
    heading = 'reaction pH',
    name='Reaction pH',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'
    )

class Command(BaseCommand):
    help='Ports database from pre-0.02 to 0.02'

    def handle(self, folder, *args, **kwargs):
        with open(path.join(folder, 'User.tsv')) as userFile:
            reader = csv.DictReader(userFile, delimiter='\t')
            for r in reader:
                u = User(
                        username=r['username'],
                        first_name=r['first_name'],
                        last_name=r['last_name'],
                        email=r['email'],
                        is_staff=int(r['is_staff']),
                        is_superuser=int(r['is_superuser']),
                    )
                u.password=r['password']
                u.save()
        with open(path.join(folder, 'labGroup.tsv')) as labGroups:
            reader = csv.DictReader(labGroups, delimiter='\t')
            for r in reader:
                l = LabGroup(**r)
                l.save()
        with open(path.join(folder, 'labgroup_users.tsv')) as labGroupUsers:
            reader = csv.DictReader(labGroupUsers, delimiter='\t')
            for r in reader:
                l = LabGroup.objects.get(title=r['title'])
                l.users.add(User.objects.get(username=r['username']))
        with open(path.join(folder, 'compound_labs.tsv')) as compounds:
            reader = csv.DictReader(compounds, delimiter='\t')
            cs = ChemSpider(settings.CHEMSPIDER_TOKEN)
            for r in reader:
                l = LabGroup.objects.get(title=r['labGroup.title'])
                CASResults = cs.simple_search(r['CAS_ID']
                nameResults = cs.simple_search(r['name'])
                if len(CASResults) != 1:
                    if len(nameResults) != 1:
                        raise RunTimeError('Could not get unambiguous chemspider entry for CAS ID {} with name{}'.format(r['CAS_ID'], r['name']))
                    else:
                        c = Compound(CSID=nameResults[0].csid, labGroup=l)
                else:
                    c = Compound(CSID=nameResults[0].csid, labGroup=l)
                c.csConsistencyCheck()
                c.save()
        with open(path.join(folder, 'chemicalClasses.tsv')) as chemicalClasses:
            reader = csv.DictReader(chemicalClasses, delimiter='\t')
            for r in reader:
                c1ChemicalClass.get_or_create(r['chemicalClass.label'])
                c2 = Compound.get(abbrev=r['compound.abbrev'])
                c2.chemicalClasses.add(c1)
        with open(path.join(folder, 'performedReactions.tsv')) as reactions:
            reader = csv.DictReader(reactions, delimiter='\t')
            for r in reader:
                p = PerformedReaction(
                    reference = r['reference'],
                    labGroup = LabGroup.get(title=r['labGroup.title']),
                    notes = r['notes'],
                    user = User.get(username=r['user.username']),
                    valid = int(r['valid']),
                    legacyRecommendedFlag=r['legacyRecommendedFlag']=='Yes',
                    insertedDateTime=r['insertedDateTime'],
                    public=int(r['public']),
                    )
                p.save()
        with open(path.join(folder, 'performedReactions.tsv')) as reactions:
            reader = csv.DictReader(reactions, delimiter='\t')
            for r in reader:
                p = PerformedReaction.get(reference=r['reference'])
                p.duplicateOf = PerformedReaction.get(reference=r['duplicate_of'])
                p.save()
                outValue = OrdRxnDescriptorValue.objects.get_or_create(descriptor=outcomeDescriptor, reaction=p)
                outValue.value = int(r['outcome']) if r['outcome'] in (str(x) for x in range (1, 5)) else None 
                outValue.save()
                outBoolValue = BoolRxnDescriptorValue.objects.get_or_create(descriptor=outcomeBooleanDescriptor, reaction=p)
                outBoolValue.value = True if outValue.value > 2 else False
                outBoolValue.save()
                purityValue = OrdRxnDescriptorValue.objects.get_or_create(descriptor=purityDescriptor, reaction=p)
                purityValue.value = int(r['purity']) if purity in ('1', '2') else None
                purityValue.save()
                temperatureDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=temperatureDescriptor, reaction=p)
                temperatureDescriptorValue.value = float(r['temp']) + 273.15
                temperatureDescriptorValue.save()
                timeDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=timeDescriptor, reaction=p)
                timeDescriptorValue.value = float(r['time'])*60
                timeDescriptorValue.save()
                pHDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=pHDescriptor, reaction=p)
                pHDescriptorValue.value = float(r['pH'])
                pHDescriptorValue.save()
        with open(path.join(folder, 'compoundquantities.tsv')) as cqs:
            reader = csv.DictReader(cqs, delimiter='\t')
            for r in reader:
                compound = Compound.objects.get(abbrev=r['compound.abbrev'])
                mw = NumMolDescriptorValue.objects.get(compound=compound, descriptor__heading='mw').value
                reaction = PerformedReaction.objects.get(reference=r['reaction.reference'])
                compoundrole = CompoundRole.objects.get_or_create(label=r['compoundrole.name'])
                if r['unit'] == 'g':
                    amount = float(r['amount'])/mw
                elif r['unit'] == 'd':
                    amount = float(r['amount'])*0.0375*float(r['density'])/mw
                elif r['unit'] == 'mL':
                    amount = float(r['amount'])*0.0375*float(r['density'])/mw
                else:
                    raise RuntimeError('invalid unit entered')
                amount = amount * 1000
                quantity = CompoundQuantity(amount=amount, role=compoundrole, compound=compound, reaction=reaction)
                quantity.save()
