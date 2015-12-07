from django.core.management.base import BaseCommand
from django.core.exceptions import ValidationError
from django.db import transaction
from django.conf import settings
from os import path
import csv
from django.contrib.auth.models import User
from DRP.models import LabGroup, Compound, PerformedReaction 
from DRP.models import ChemicalClass, NumRxnDescriptor
from DRP.models import BoolRxnDescriptor, OrdRxnDescriptor
from DRP.models import NumRxnDescriptorValue, BoolRxnDescriptorValue
from DRP.models import OrdRxnDescriptorValue, NumMolDescriptorValue 
from DRP.models import CompoundRole, CompoundQuantity
from chemspipy import ChemSpider

outcomeDescriptor = OrdRxnDescriptor.objects.get_or_create(
    heading='crystallisation_outcome',
    name='Four class crystallisation outcome',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0',
    maximum=4,
    minimum=1
    )[0]

outcomeBooleanDescriptor = BoolRxnDescriptor.objects.get_or_create(
    heading='boolean_crystallisation_outcome',
    name='Two class crystallisation outcome',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'    
    )[0]

purityDescriptor = OrdRxnDescriptor.objects.get_or_create(
    heading='crystallisation_purity_outcome',
    name='Two class crystallisation purity',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0',
    maximum=2,
    minimum=1
    )[0]

temperatureDescriptor = NumRxnDescriptor.objects.get_or_create(
    heading='reaction_temperature',
    name='Reaction temperature in degrees C',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0' 
    )[0]

timeDescriptor = NumRxnDescriptor.objects.get_or_create(
    heading='reaction_time',
    name='Reaction time in minutes',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'  
    )[0]

pHDescriptor = NumRxnDescriptor.objects.get_or_create(
    heading = 'reaction_pH',
    name='Reaction pH',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'
    )[0]

preHeatStandingDescriptor = NumRxnDescriptor.objects.get_or_create(
    heading='pre_heat_standing',
    name='Pre heat standing time',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'
    )[0]

teflonDescriptor = BoolRxnDescriptor.objects.get_or_create(
    heading='teflon_pouch',
    name='Was this reaction performed in a teflon pouch?',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'
    )[0]

class Command(BaseCommand):
    help='Ports database from pre-0.02 to 0.02'

    def handle(self, folder, *args, **kwargs):
        with open(path.join(folder, 'User.tsv')) as userFile:
            reader = csv.DictReader(userFile, delimiter='\t')
            for r in reader:
                if not User.objects.filter(username=r['username']).exists():
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
                if not LabGroup.objects.filter(title=r['title']).exists():
                    l = LabGroup(**r)
                    l.save()
        with open(path.join(folder, 'labgroup_users.tsv')) as labGroupUsers:
            reader = csv.DictReader(labGroupUsers, delimiter='\t')
            for r in reader:
                l = LabGroup.objects.get(title=r['title'])
                l.users.add(User.objects.get(username=r['username']))
        with open(path.join(folder, 'performedReactions.tsv')) as reactions:
            reader = csv.DictReader(reactions, delimiter='\t')
            for r in reader:
                if not PerformedReaction.objects.filter(reference=r['reference'].lower()).exists():
                    p = PerformedReaction(
                        reference = r['reference'],
                        labGroup = LabGroup.objects.get(title=r['labGroup.title']),
                        notes = r['notes'],
                        user = User.objects.get(username=r['user.username']),
                        valid = int(r['valid']),
                        legacyRecommendedFlag=r['legacyRecommendedFlag']=='Yes',
                        insertedDateTime=r['insertedDateTime'],
                        public=int(r['public'])
                        )
                    self.stdout.write('Creating reaction with reference {}'.format(p.reference))
                    p.validate_unique()
                    p.save(calcDescriptors=False)
        with open(path.join(folder, 'performedReactions.tsv')) as reactions:
            reader = csv.DictReader(reactions, delimiter='\t')
            for r in reader:
                self.stdout.write('Reiterating for reaction with reference {}'.format(r['reference'].lower()))
                ps = PerformedReaction.objects.filter(reference=r['reference'].lower())
                if ps.count() > 1:
                    ps = ps.filter(valid=True)
                if ps.exists():
                    if ps.count() > 1:
                        raise RuntimeError('{} has more than one reaction'.format(r['reference'].lower()))
                    p = ps[0]
                    try:
                        p.duplicateOf = PerformedReaction.objects.get(reference=r['duplicateOf.reference'].lower())
                        p.save()
                    except PerformedReaction.DoesNotExist:
                        pass
                    outValue = OrdRxnDescriptorValue.objects.get_or_create(descriptor=outcomeDescriptor, reaction=p)[0]
                    outValue.value = int(r['outcome']) if (r['outcome'] in (str(x) for x in range (1, 5))) else None 
                    outValue.save()
                    outBoolValue = BoolRxnDescriptorValue.objects.get_or_create(descriptor=outcomeBooleanDescriptor, reaction=p)[0]
                    outBoolValue.value = True if (outValue.value > 2) else False
                    outBoolValue.save()
                    purityValue = OrdRxnDescriptorValue.objects.get_or_create(descriptor=purityDescriptor, reaction=p)[0]
                    purityValue.value = int(r['purity']) if (r['purity'] in ('1', '2')) else None
                    purityValue.save()
                    temperatureDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=temperatureDescriptor, reaction=p)[0]
                    temperatureDescriptorValue.value = (float(r['temp']) + 273.15) if (r['temp'] not in ('', '?')) else None
                    temperatureDescriptorValue.save()
                    timeDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=timeDescriptor, reaction=p)[0]
                    timeDescriptorValue.value = float(r['time'])*60 if (r['time'] not in ['', '?']) else None
                    timeDescriptorValue.save()
                    pHDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=pHDescriptor, reaction=p)[0]
                    pHDescriptorValue.value = float(r['pH']) if (r['pH'] not in ('', '?')) else None
                    pHDescriptorValue.save()
                    preHeatStandingDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=preHeatStandingDescriptor, reaction=p)[0]
                    preHeatStandingDescriptorValue.value = bool(r['pre_heat standing']) if (r.get('pre_heat standing') not in ('', None)) else None
                    preHeatStandingDescriptorValue.save()
                    teflonDescriptorValue = BoolRxnDescriptorValue.objects.get_or_create(descriptor=teflonDescriptor, reaction=p)[0]
                    teflonDescriptorValue.value = bool(int(r['teflon_pouch'])) if (r.get('teflon_pouch') not in(None, '')) else None
                    teflonDescriptorValue.save()
        with open(path.join(folder, 'compound_labs.tsv')) as compounds:
            reader = csv.DictReader(compounds, delimiter='\t')
            cs = ChemSpider(settings.CHEMSPIDER_TOKEN)
            for r in reader:
                l = LabGroup.objects.get(title=r['labGroup.title'])
                if not Compound.objects.filter(abbrev=r['abbrev']).exists():
                    self.stdout.write('Importing compound with abbreviation {} and name {}'.format(r['abbrev'], r['name']))
                    if r.get('custom') != '1':
                        try: 
                            if r.get('CSID') not in ('', None):
                                c = Compound(CSID=r['CSID'], labGroup=l, abbrev=r['abbrev'])
                                c.csConsistencyCheck()
                                c.save()
                            else:
                                if r.get('CAS_ID') not in (None, ''):
                                    CASResults = cs.simple_search(r['CAS_ID'])
                                else:
                                    CASResults = []
                                if len(CASResults) != 1:
                                        nameResults = cs.simple_search(r.get('name'))
                                        if len(nameResults) != 1:
                                            raise RuntimeError('Could not get unambiguous chemspider entry for CAS ID {} with name{}'.format(r['CAS_ID'], r['name']))
                                        else:
                                            c = Compound(CSID=nameResults[0].csid, labGroup=l, abbrev=r['abbrev'])
                                            c.csConsistencyCheck()
                                            c.save()
                                else:
                                    c = Compound(CSID=CASResults[0].csid, labGroup=l, abbrev=r['abbrev'])
                                    c.csConsistencyCheck()
                                    c.save()
                        except ValidationError as e:
                            c.delete()
                            raise e
                    else:
                        if r.get('INCHI') is None:
                            r['INCHI'] = ''
                        if r.get('smiles') is None:
                            r['smiles'] = ''
                        c = Compound.objects.get_or_create(labGroup=l, custom=True, name=r['name'], abbrev=r['abbrev'], formula=r['formula'], smiles=r['smiles'], INCHI=r['INCHI'])[0]
                    self.stdout.write(c.name.encode('utf-8'))
                    c.save()
        with open(path.join(folder, 'compound_chemicalClasses.tsv')) as chemicalClasses:
            reader = csv.DictReader(chemicalClasses, delimiter='\t')
            for r in reader:
                self.stdout.write('working with class {}'.format(r['chemicalClass.label']))
                cs = Compound.objects.filter(abbrev=r['compound.abbrev'])
                if cs.count() > 0:
                    c1=ChemicalClass.objects.get_or_create(label=r['chemicalClass.label'])[0]
                    for c2 in cs:
                        if not c1 in c2.chemicalClasses.all():
                            c2.chemicalClasses.add(c1)
                            c2.save()
        with open(path.join(folder, 'compoundquantities.tsv')) as cqs:
            reader = csv.DictReader(cqs, delimiter='\t')
            for r in reader:
                try:
                    reaction = PerformedReaction.objects.get(reference=r['reaction.reference'])
                    compound = Compound.objects.get(abbrev=r['compound.abbrev'], labGroup=reaction.labGroup)
                    if r['compound.abbrev'] in ('water', 'H2O'):
                        r['density'] = 1
                    mw = NumMolDescriptorValue.objects.get(compound=compound, descriptor__heading='mw').value
                    if r['compoundrole.name'] != 'pH':
                        self.stdout.write('adding {} to {}'.format(compound.abbrev, reaction.reference))
                        compoundrole = CompoundRole.objects.get_or_create(label=r['compoundrole.name'])[0]
                        if r['amount'] in ('', '?'):
                            amount = None
                        elif r['unit'] == 'g':
                            amount = float(r['amount'])/mw
                        elif r['unit'] == 'd':
                            amount = float(r['amount'])*0.0375*float(r['density'])/mw
                        elif r['unit'] == 'mL':
                            amount = float(r['amount'])*float(r['density'])/mw
                        else:
                            raise RuntimeError('invalid unit entered')
                        if amount is not None:
                            amount = (amount * 1000)
                        cqq = CompoundQuantity.objects.filter(role=compoundrole, compound=compound, reaction=reaction)
                        if cqq.count() > 1:
                            cqq.delete()
                        quantity = CompoundQuantity.objects.get_or_create(role=compoundrole, compound=compound, reaction=reaction)[0]
                        quantity.amount = amount
                        quantity.save()
                    else:
                        reaction.notes += ' pH adjusting reagent used: {}, {}{}'.format(r['compound.abbrev'], r['amount'], r['unit'])
                        reaction.save(calcDescriptors=False)
                except Compound.DoesNotExist as e:
                    reaction.notes += ' Unknown Reactant {} with amount {} {}'.format(r['compound.abbrev'], r['amount'], r['unit'])
                    reaction.valid = False
                    reaction.save(calcDescriptors = False)
                except PerformedReaction.DoesNotExist as e:
                    pass
