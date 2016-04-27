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

leakDescriptor = BoolRxnDescriptor.objects.get_or_create(
    heading='leak',
    name='Did this reaction leak?',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'
    )[0]

slowCoolDescriptor = BoolRxnDescriptor.objects.get_or_create(
    heading='slow_cool',
    name='Was a slow cool performed for this reaction?',
    calculatorSoftware='manual',
    calculatorSoftwareVersion='0'
    )[0]

class Command(BaseCommand):
    help='Ports database from pre-0.02 to 0.02'

    def add_arguments(self, parser):
        parser.add_argument('directory', help='The directory where the tsv files are')

    def handle(self, *args, **kwargs):
        folder = kwargs['directory']
        
        with open(path.join(folder, 'performedReactions.tsv')) as reactions:
            reader = csv.DictReader(reactions, delimiter='\t')
            for r in reader:
                ref = r['reference'].lower()
                ps = PerformedReaction.objects.filter(reference=ref)
                if ps.exists():
                    if ps.count() > 1:
                        raise RuntimeError('More than one reaction with reference {}'.fromat(ref))
                    self.stdout.write('Updating reaction with reference {}'.format(ref))
                    ps.update(labGroup = LabGroup.objects.get(title=r['labGroup.title']),
                                notes = r['notes'],
                                user = User.objects.get(username=r['user.username']),
                                valid = int(r['valid']),
                                legacyRecommendedFlag=r['legacyRecommendedFlag']=='Yes',
                                insertedDateTime=r['insertedDateTime'],
                                public=int(r['public'])
                              )
                else:
                    p = PerformedReaction(
                        reference = ref,
                        labGroup = LabGroup.objects.get(title=r['labGroup.title']),
                        notes = r['notes'],
                        user = User.objects.get(username=r['user.username']),
                        valid = int(r['valid']),
                        legacyRecommendedFlag=r['legacyRecommendedFlag']=='Yes',
                        insertedDateTime=r['insertedDateTime'],
                        public=int(r['public'])
                        )
                    self.stdout.write('Creating reaction with reference {}'.format(ref))
                    p.validate_unique()
                    p.save(calcDescriptors=False)
        with open(path.join(folder, 'performedReactions.tsv')) as reactions:
            reader = csv.DictReader(reactions, delimiter='\t')
            outValues = []
            outBoolValues = []
            purityValues = []
            temperatureValues = []
            timeValues = []
            pHValues = []
            preHeatStandingValues = []
            teflonValues = []
            slowCoolValues = []
            leakValues = []
            
            for r in reader:
                ref = r['reference'].lower()
                self.stdout.write('Reiterating for reaction with reference {}'.format(ref))
                ps = PerformedReaction.objects.filter(reference=ref)
                if not ps:
                    raise RuntimeError('Cannot find reaction with reference {}'.format(ref)
                else:
                    if ps.count() > 1:
                        raise RuntimeError('More than one reaction with reference {}'.format(ref))
                    p = ps[0]
                    try:
                        p.duplicateOf = PerformedReaction.objects.get(reference=r['duplicateOf.reference'].lower())
                        p.save()
                    except PerformedReaction.DoesNotExist:
                        pass

                    outcomeValue = int(r['outcome']) if (r['outcome'] in (str(x) for x in range (1, 5))) else None
                    try:
                        v = OrdRxnDescriptorValue.objects.get(descriptor=outcomeDescriptor, reaction=p)
                        if v.value != outcomeValue:
                            v.value = outcomeValue
                            v.save()
                    except OrdRxnDescriptorValue.DoesNotExist:
                        outValue = outcomeDescriptor.createValue(p, outcomeValue)
                        outValues.append(outValue)
                    
                    value = True if (outcomeValue > 2) else False
                    try:
                        v = BoolRxnDescriptorValue.objects.get(descriptor=outcomeBooleanDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except BoolRxnDescriptorValue.DoesNotExist:
                        outBoolValue = outcomeBooleanDescriptor.createValue(p, value)
                        outBoolValues.append(outBoolValue)
                    
                    value = int(r['purity']) if (r['purity'] in ('1', '2')) else None
                    try:
                        v = OrdRxnDescriptorValue.objects.get(descriptor=purityDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except OrdRxnDescriptorValue.DoesNotExist:
                        purityValue = purityDescriptor.createValue(p, value)
                        purityValues.append(purityValue)
                    
                    value = (float(r['temp']) + 273.15) if (r['temp'] not in ('', '?')) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=temperatureDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        temperatureDescriptorValue = temperatureDescriptor.createValue(p, value)
                        temperatureValues.append(temperatureDescriptorValue)
                    
                    value = float(r['time'])*60 if (r['time'] not in ['', '?']) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=timeDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        timeDescriptorValue = timeDescriptor.createValue(p, value)
                        timeValues.append(timeDescriptorValue)
                    
                    value = float(r['pH']) if (r['pH'] not in ('', '?')) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=pHDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        pHDescriptorValue = pHDescriptor.createValue(p, value)
                        pHValues.append(pHDescriptorValue)
                    
                    value = bool(r['pre_heat standing']) if (r.get('pre_heat standing') not in ('', None)) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=preHeatStandingDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        preHeatStandingDescriptorValue = preHeatStandingDescriptor.createValue(p, value)
                        preHeatStandingValues.append(preHeatStandingDescriptorValue)
                    
                    value = bool(int(r['teflon_pouch'])) if (r.get('teflon_pouch') not in(None, '')) else None
                    try:
                        v = BoolRxnDescriptorValue.objects.get(descriptor=teflonDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except BoolRxnDescriptorValue.DoesNotExist:
                        teflonDescriptorValue = teflonDescriptor.createValue(p, value)
                        teflonValues.append(teflonDescriptorValue)

                    leak_string = r['leak']
                    if leak_string in (None, ''):
                        value = None
                    elif leak_string.lower() == 'yes':
                        value = True
                    elif leak_string.lower() == 'no':
                        value = False
                    else:
                        raise RuntimeError("Unrecognized string '{}' in leak column".format(leak_string))
                    try:
                        v = BoolRxnDescriptorValue.objects.get(descriptor=leakDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except BoolRxnDescriptorValue.DoesNotExist:
                        leakDescriptorValue = teflonDescriptor.createValue(p, value)
                        leakValues.append(leakDescriptorValue)

                    slow_cool_string = r['slow_cool']
                    if slow_cool_string in (None, ''):
                        value = None
                    elif slow_cool_string.lower() == 'yes':
                        value = True
                    elif slow_cool_string.lower() == 'no':
                        value = False
                    else:
                        raise RuntimeError("Unrecognized string '{}' in slow_cool column".format(slow_cool_string))
                    try:
                        v = BoolRxnDescriptorValue.objects.get(descriptor=slowCoolDescriptor, reaction=p)
                        if v.value != value:
                            v.value = value
                            v.save()
                    except BoolRxnDescriptorValue.DoesNotExist:
                        slowCoolDescriptorValue = slowCoolDescriptor.createValue(p, value)
                        slowCoolValues.append(slowCoolDescriptorValue)

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
                        BoolRxnDescriptorValue.objects.bulk_create(leakValues)
                        BoolRxnDescriptorValue.objects.bulk_create(slowCoolValues)
                        
                        
                        outValues = []
                        outBoolValues = []
                        purityValues = []
                        temperatureValues = []
                        timeValues = []
                        pHValues = []
                        preHeatStandingValues = []
                        teflonValues = []
                        leakValues = []
                        slowCoolValues = []
                        self.stdout.write("...saved")
                    
        with open(path.join(folder, 'compoundquantities.tsv')) as cqs:
            reader = csv.DictReader(cqs, delimiter='\t')
            quantities = []
            for r in reader:
                try:
                    reaction = PerformedReaction.objects.get(reference=r['reaction.reference'].lower())
                    compound = Compound.objects.get(abbrev=r['compound.abbrev'], labGroup=reaction.labGroup)
                    if r['compound.abbrev'] in ('water', 'H2O'):
                        r['density'] = 1
                    mw = NumMolDescriptorValue.objects.get(compound=compound, descriptor__heading='mw').value
                    if r['compoundrole.name'] in (None, '', '?'):
                        self.stderr.write('No role for reactant {} with amount {} {} in reaction {}'.format(r['compound.abbrev'], r['amount'], r['unit'], r['reaction.reference']))
                        reaction.notes += ' No role for reactant {} with amount {} {}'.format(r['compound.abbrev'], r['amount'], r['unit'])
                        reaction.save(calcDescriptors=False)
                    elif r['compoundrole.name'] == 'pH':
                        reaction.notes += ' pH adjusting reagent used: {}, {}{}'.format(r['compound.abbrev'], r['amount'], r['unit'])
                        reaction.save(calcDescriptors=False)
                    else:
                        self.stdout.write('adding {} to {}'.format(compound.abbrev, reaction.reference))
                        compoundrole = CompoundRole.objects.get_or_create(label=r['compoundrole.name'])[0]
                        if r['amount'] in ('', '?'):
                            amount = None
                            reaction.notes += ' No amount for reactant {} with role {}'.format(r['compound.abbrev'], r['compound.role'])
                        elif r['unit'] == 'g':
                            amount = float(r['amount'])/mw
                        elif r['unit'] == 'd':
                            amount = float(r['amount'])*0.0375*float(r['density'])/mw
                        elif r['unit'] == 'mL':
                            amount = float(r['amount'])*float(r['density'])/mw
                        else:
                            raise RuntimeError('invalid unit entered')
                        # convert to milimoles
                        if amount is not None:
                            amount = (amount * 1000)
                        cqq = CompoundQuantity.objects.filter(compound=compound, reaction=reaction)
                        if cqq.count() > 1:
                            cqq.delete()
                        quantity = CompoundQuantity.objects.get_or_create(role=compoundrole, compound=compound, reaction=reaction)[0]
                        quantity.amount = amount
                        quantity.save()

                        try:
                            quantity = CompoundQuantity.objects.get(compound=compound, reaction=reaction)
                            if quantity.amount != amount or quantity.role != compoundrole:
                                quantity.amount = amount
                                quantity.role = compoundrole
                        except CompoundQuantity.DoesNotExist:
                            quantity = CompoundQuantity(compound=compound, reaction=reaction, role=compoundrole)
                            quantities.append(quantity)

                        if len(quantities) > 500:
                            CompoundQuantity.objects.bulk_create(quantities)
                except Compound.DoesNotExist as e:
                    self.stderr.write('Unknown Reactant {} with amount {} {} in reaction {}'.format(r['compound.abbrev'], r['amount'], r['unit'], r['reaction.reference']))
                    reaction.notes += ' Unknown Reactant {} with amount {} {}'.format(r['compound.abbrev'], r['amount'], r['unit'])
                    reaction.valid = False
                    reaction.save(calcDescriptors=False)
                except PerformedReaction.DoesNotExist as e:
                    self.stderr.write('Cannot find reactions {}'.format(r['reaction.reference']))
                    raise e
