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
    heading='reaction_pH',
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
    help = 'Ports database from pre-0.02 to 0.02'

    def add_arguments(self, parser):
        parser.add_argument('directory', help='The directory where the tsv files are')

    def handle(self, *args, **kwargs):
        folder = kwargs['directory']
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
                    u.password = r['password']
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

        if not path.isfile(path.join(folder, 'performedReactionsNoDups.tsv')):
            self.stdout.write('Writing file with duplicate references disambiguated (arbitrarily)')
            with open(path.join(folder, 'performedReactions.tsv')) as in_file, open(path.join(folder, 'performedReactionsNoDups.tsv'), 'w') as out_file:
                references = set()
                reader = csv.DictReader(in_file, delimiter='\t')
                writer = csv.DictWriter(out_file, delimiter='\t', fieldnames=reader.fieldnames)
                writer.writeheader()

                case_count = 0
                valid_case_count = 0
                dup_count = 0
                for r in reader:
                    ref = r['reference'].lower()
                    if ref != r['reference']:
                        self.stderr.write('Reference {} was not in lowercase. Converted.'.format(r['reference']))
                        case_count += 1
                        if r['valid'] == '1':
                            valid_case_count += 1

                    if ref in references:
                        r['notes'] += ' Duplicated reference'
                        r['valid'] = 0
                        dup_count += 1
                        i = 1
                        new_ref = ref
                        while new_ref in references:
                            new_ref = '{}_dup{}'.format(ref, i)
                            i += 1
                        self.stderr.write('Reference {} duplicated {} times. Renamed and invalidated'.format(ref, i))
                        ref = new_ref
                    references.add(ref)
                    r['reference'] = ref
                    writer.writerow(r)
            self.stderr.write('{} references converted to lowercase. {} were valid'.format(case_count, valid_case_count))
            self.stderr.write('{} references with _dupX appended to remove duplicate reference'.format(dup_count))

        with open(path.join(folder, 'performedReactionsNoDups.tsv')) as reactions:
            reader = csv.DictReader(reactions, delimiter='\t')
            for r in reader:
                if not PerformedReaction.objects.filter(reference=r['reference'].lower()).exists():
                    p = PerformedReaction(
                        reference=r['reference'],
                        labGroup=LabGroup.objects.get(title=r['labGroup.title']),
                        notes=r['notes'],
                        user=User.objects.get(username=r['user.username']),
                        valid=int(r['valid']),
                        legacyRecommendedFlag=r['legacyRecommendedFlag'] == 'Yes',
                        insertedDateTime=r['insertedDateTime'],
                        public=int(r['public'])
                    )
                    self.stdout.write('Creating reaction with reference {}'.format(p.reference))
                    p.validate_unique()
                    p.save(calcDescriptors=False)
        with open(path.join(folder, 'performedReactionsNoDups.tsv')) as reactions:
            reader = csv.DictReader(reactions, delimiter='\t')
            outValues = []
            outBoolValues = []
            purityValues = []
            temperatureValues = []
            timeValues = []
            pHValues = []
            preHeatStandingValues = []
            teflonValues = []

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

                    #outValue = OrdRxnDescriptorValue.objects.get_or_create(descriptor=outcomeDescriptor, reaction=p)[0]
                    outcomeValue = int(r['outcome']) if (r['outcome'] in (str(x) for x in range(1, 5))) else None
                    try:
                        v = OrdRxnDescriptorValue.objects.get(descriptor=outcomeDescriptor, reaction=p)
                        if v.value != outcomeValue:
                            v.value = outcomeValue
                            v.save()
                    except OrdRxnDescriptorValue.DoesNotExist:
                        outValue = outcomeDescriptor.createValue(p, outcomeValue)
                        # outValue.save()
                        outValues.append(outValue)

                    #outBoolValue = BoolRxnDescriptorValue.objects.get_or_create(descriptor=outcomeBooleanDescriptor, reaction=p)[0]
                    value = True if (outcomeValue > 2) else False
                    try:
                        v = BoolRxnDescriptorValue.objects.get(descriptor=outcomeBooleanDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except BoolRxnDescriptorValue.DoesNotExist:
                        # outBoolValue.save()
                        outBoolValue = outcomeBooleanDescriptor.createValue(p, value)
                        outBoolValues.append(outBoolValue)

                    #purityValue = OrdRxnDescriptorValue.objects.get_or_create(descriptor=purityDescriptor, reaction=p)[0]
                    value = int(r['purity']) if (r['purity'] in ('1', '2')) else None
                    try:
                        v = OrdRxnDescriptorValue.objects.get(descriptor=purityDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except OrdRxnDescriptorValue.DoesNotExist:
                        # purityValue.save()
                        purityValue = purityDescriptor.createValue(p, value)
                        purityValues.append(purityValue)

                    #temperatureDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=temperatureDescriptor, reaction=p)[0]
                    value = (float(r['temp']) + 273.15) if (r['temp'] not in ('', '?')) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=temperatureDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        # temperatureDescriptorValue.save()
                        temperatureDescriptorValue = temperatureDescriptor.createValue(p, value)
                        temperatureValues.append(temperatureDescriptorValue)

                    #timeDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=timeDescriptor, reaction=p)[0]
                    value = float(r['time']) * 60 if (r['time'] not in ['', '?']) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=timeDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        # timeDescriptorValue.save()
                        timeDescriptorValue = timeDescriptor.createValue(p, value)
                        timeValues.append(timeDescriptorValue)

                    #pHDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=pHDescriptor, reaction=p)[0]
                    value = float(r['pH']) if (r['pH'] not in ('', '?')) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=pHDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        # pHDescriptorValue.save()
                        pHDescriptorValue = pHDescriptor.createValue(p, value)
                        pHValues.append(pHDescriptorValue)

                    #preHeatStandingDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=preHeatStandingDescriptor, reaction=p)[0]
                    value = bool(r['pre_heat standing']) if (r.get('pre_heat standing') not in ('', None)) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=preHeatStandingDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        # preHeatStandingDescriptorValue.save()
                        preHeatStandingDescriptorValue = preHeatStandingDescriptor.createValue(p, value)
                        preHeatStandingValues.append(preHeatStandingDescriptorValue)

                    #teflonDescriptorValue = BoolRxnDescriptorValue.objects.get_or_create(descriptor=teflonDescriptor, reaction=p)[0]
                    value = bool(int(r['teflon_pouch'])) if (r.get('teflon_pouch') not in(None, '')) else None
                    try:
                        v = BoolRxnDescriptorValue.objects.get(descriptor=teflonDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except BoolRxnDescriptorValue.DoesNotExist:
                        # teflonDescriptorValue.save()
                        teflonDescriptorValue = teflonDescriptor.createValue(p, value)
                        teflonValues.append(teflonDescriptorValue)

                    if len(outValues) > 500:
                        self.stdout.write("Saving...")
                        OrdRxnDescriptorValue.objects.bulk_create(outValues)
                        BoolRxnDescriptorValue.objects.bulk_create(outBoolValues)
                        OrdRxnDescriptorValue.objects.bulk_create(purityValues)
                        NumRxnDescriptorValue.objects.bulk_create(temperatureValues)
                        NumRxnDescriptorValue.objects.bulk_create(timeValues)
                        NumRxnDescriptorValue.objects.bulk_create(pHValues)
                        NumRxnDescriptorValue.objects.bulk_create(preHeatStandingValues)
                        BoolRxnDescriptorValue.objects.bulk_create(teflonValues)

                        outValues = []
                        outBoolValues = []
                        purityValues = []
                        temperatureValues = []
                        timeValues = []
                        pHValues = []
                        preHeatStandingValues = []
                        teflonValues = []
                        self.stdout.write("...saved")

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
                    c1 = ChemicalClass.objects.get_or_create(label=r['chemicalClass.label'])[0]
                    for c2 in cs:
                        if not c1 in c2.chemicalClasses.all():
                            c2.chemicalClasses.add(c1)
                            c2.save()
        with open(path.join(folder, 'compoundquantities.tsv')) as cqs:
            reader = csv.DictReader(cqs, delimiter='\t')
            for r in reader:
                try:
                    reaction = PerformedReaction.objects.get(reference=r['reaction.reference'].lower())
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
                            amount = float(r['amount']) / mw
                        elif r['unit'] == 'd':
                            amount = float(r['amount']) * 0.0375 * float(r['density']) / mw
                        elif r['unit'] == 'mL':
                            amount = float(r['amount']) * float(r['density']) / mw
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
                    self.stderr.write('Unknown Reactant {} with amount {} {} in reaction {}'.format(r['compound.abbrev'], r['amount'], r['unit'], r['reaction.reference']))
                    raw_input("Continue?")
                    reaction.notes += ' Unknown Reactant {} with amount {} {}'.format(r['compound.abbrev'], r['amount'], r['unit'])
                    reaction.valid = False
                    reaction.save(calcDescriptors=False)
                except PerformedReaction.DoesNotExist as e:
                    raise e
