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
import warnings

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

    def add_arguments(self, parser):
        parser.add_argument('directory', help='The directory where the tsv files are')

    def handle(self, *args, **kwargs):
        folder = kwargs['directory']
        if not path.isfile(path.join(folder, 'performedReactionsNoDupsLower.tsv')):
            self.stdout.write('Writing file with all references that were uppercase (now lower) and duplicate references disambiguated (arbitrarily)')
            with open(path.join(folder, 'performedReactions.tsv')) as in_file, open(path.join(folder, 'performedReactionsNoDupsLower.tsv'), 'w') as out_file:
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
        with open(path.join(folder, 'performedReactionsNoDupsLower.tsv')) as reactions:
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
                ref = r['reference'].lower()
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
                    outcomeValue = int(r['outcome']) if (r['outcome'] in (str(x) for x in range (1, 5))) else None
                    try:
                        v = OrdRxnDescriptorValue.objects.get(descriptor=outcomeDescriptor, reaction=p)
                        if v.value != outcomeValue:
                            sys.stderr.write('Value for reaction {} and descriptor {} found with different value. Overwriting old value {} with new value {}'.format(ref, outcomeDescriptor, v.value, outcomeValue))
                            v.value = outcomeValue
                            v.save()
                    except OrdRxnDescriptorValue.DoesNotExist:
                        outValue = outcomeDescriptor.createValue(p, outcomeValue)
                        #outValue.save()
                        outValues.append(outValue)
                    
                    #outBoolValue = BoolRxnDescriptorValue.objects.get_or_create(descriptor=outcomeBooleanDescriptor, reaction=p)[0]
                    value = True if (outcomeValue > 2) else False
                    try:
                        v = BoolRxnDescriptorValue.objects.get(descriptor=outcomeBooleanDescriptor, reaction=p)
                        if v.value != value:
                            sys.stderr.write('Value for reaction {} and descriptor {} found with different value. Overwriting old value {} with new value {}'.format(ref, outcomeBooleanDescriptor, v.value, value))
                            v.value = value
                            v.save()
                    except BoolRxnDescriptorValue.DoesNotExist:
                        #outBoolValue.save()
                        outBoolValue = outcomeBooleanDescriptor.createValue(p, value)
                        outBoolValues.append(outBoolValue)
                    
                    #purityValue = OrdRxnDescriptorValue.objects.get_or_create(descriptor=purityDescriptor, reaction=p)[0]
                    value = int(r['purity']) if (r['purity'] in ('1', '2')) else None
                    try:
                        v = OrdRxnDescriptorValue.objects.get(descriptor=purityDescriptor, reaction=p)
                        if v.value != value:
                            sys.stderr.write('Value for reaction {} and descriptor {} found with different value. Overwriting old value {} with new value {}'.format(ref, purityDescriptor, v.value, value))
                            v.value = value
                            v.save()
                    except OrdRxnDescriptorValue.DoesNotExist:
                        #purityValue.save()
                        purityValue = purityDescriptor.createValue(p, value)
                        purityValues.append(purityValue)
                    
                    #temperatureDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=temperatureDescriptor, reaction=p)[0]
                    value = (float(r['temp']) + 273.15) if (r['temp'] not in ('', '?')) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=temperatureDescriptor, reaction=p)
                        if v.value != value:
                            sys.stderr.write('Value for reaction {} and descriptor {} found with different value. Overwriting old value {} with new value {}'.format(ref, temperatureDescriptor, v.value, value))
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        #temperatureDescriptorValue.save()
                        temperatureDescriptorValue = temperatureDescriptor.createValue(p, value)
                        temperatureValues.append(temperatureDescriptorValue)
                    
                    #timeDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=timeDescriptor, reaction=p)[0]
                    value = float(r['time'])*60 if (r['time'] not in ['', '?']) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=timeDescriptor, reaction=p)
                        if v.value != value:
                            sys.stderr.write('Value for reaction {} and descriptor {} found with different value. Overwriting old value {} with new value {}'.format(ref, timeDescriptor, v.value, value))
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        #timeDescriptorValue.save()
                        timeDescriptorValue = timeDescriptor.createValue(p, value)
                        timeValues.append(timeDescriptorValue)
                    
                    #pHDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=pHDescriptor, reaction=p)[0]
                    value = float(r['pH']) if (r['pH'] not in ('', '?')) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=pHDescriptor, reaction=p)
                        if v.value != value:
                            sys.stderr.write('Value for reaction {} and descriptor {} found with different value. Overwriting old value {} with new value {}'.format(ref, pHDescriptor, v.value, value))
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        #pHDescriptorValue.save()
                        pHDescriptorValue = pHDescriptor.createValue(p, value)
                        pHValues.append(pHDescriptorValue)
                    
                    #preHeatStandingDescriptorValue = NumRxnDescriptorValue.objects.get_or_create(descriptor=preHeatStandingDescriptor, reaction=p)[0]
                    value = bool(r['pre_heat standing']) if (r.get('pre_heat standing') not in ('', None)) else None
                    try:
                        v = NumRxnDescriptorValue.objects.get(descriptor=preHeatStandingDescriptor, reaction=p)
                        if v.value != value:
                            sys.stderr.write('Value for reaction {} and descriptor {} found with different value. Overwriting old value {} with new value {}'.format(ref, preHeatStandingDescriptor, v.value, value))
                            v.value = value
                            v.save()
                    except NumRxnDescriptorValue.DoesNotExist:
                        #preHeatStandingDescriptorValue.save()
                        preHeatStandingDescriptorValue = preHeatStandingDescriptor.createValue(p, value)
                        preHeatStandingValues.append(preHeatStandingDescriptorValue)
                    
                    #teflonDescriptorValue = BoolRxnDescriptorValue.objects.get_or_create(descriptor=teflonDescriptor, reaction=p)[0]
                    value = bool(int(r['teflon_pouch'])) if (r.get('teflon_pouch') not in(None, '')) else None
                    try:
                        v = BoolRxnDescriptorValue.objects.get(descriptor=teflonDescriptor, reaction=p)
                        if v.value != value:
                            sys.stderr.write('Value for reaction {} and descriptor {} found with different value. Overwriting old value {} with new value {}'.format(ref, teflonDescriptor, v.value, value))
                            v.value = value
                            v.save()
                    except BoolRxnDescriptorValue.DoesNotExist:
                        #teflonDescriptorValue.save()
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

                else:
                    self.stderr.write("Reaction with reference {} not found".format(r['reference'].lower()))

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
