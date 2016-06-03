import django
django.setup()
import csv
from sys import argv
from DRP.models import NumRxnDescriptor, BoolRxnDescriptor, OrdRxnDescriptor, CatRxnDescriptor, PerformedReaction, NumRxnDescriptorValue, BoolRxnDescriptorValue, OrdRxnDescriptorValue, CatRxnDescriptorValue
from DRP.management.commands import reimport_reactions

def create_descriptors(descriptor_type_dict, save=False, overwrite=False, heading_suffix='_legacy', name_prefix='Legacy descriptor with heading ', calculatorSoftware='legacy'):
    descriptor_dict = {}
    for descriptor_heading, descriptor_type in descriptor_type_dict.items():
        heading = descriptor_heading + heading_suffix
        name = name_prefix + descriptor_heading
        desc = None
        if descriptor_type == 'numeric':
            new_desc = NumRxnDescriptor(
                        heading=heading,
                        calculatorSoftware=calculatorSoftware,
                        name=name,
                        )
            try:
                desc = NumRxnDescriptor.objects.get(heading=heading)
            except NumRxnDescriptor.DoesNotExist:
                desc = new_desc
            else:
                if overwrite:
                    desc.delete()
                    desc = new_desc
                    
        elif descriptor_type == 'boolean':
            new_desc = BoolRxnDescriptor(
                        heading=heading,
                        calculatorSoftware=calculatorSoftware,
                        name=name,
                        )
            try:
                desc = BoolRxnDescriptor.objects.get(heading=heading)
            except BoolRxnDescriptor.DoesNotExist:
                desc = new_desc
            else:
                if overwrite:
                    desc.delete()
                    desc = new_desc
                    
        elif descriptor_type == 'ordinal':
            if descriptor_heading == 'purity': #hacky, but whatever
                new_desc = OrdRxnDescriptor(
                        heading=heading,
                        calculatorSoftware='legacy',
                        minimum=0,
                        maximum=2
                        )
            elif descriptor_heading == 'outcome': #hacky, but whatever
                new_desc = OrdRxnDescriptor(
                        heading=heading,
                        calculatorSoftware='legacy',
                        minimum=1,
                        maximum=4
                        )
            else:
                raise ValueError('Unrecognized ordinal descriptor')
            try:
                desc = OrdRxnDescriptor.objects.get(heading=heading)
            except OrdRxnDescriptor.DoesNotExist:
                desc = new_desc
            else:
                if overwrite:
                    desc.delete()
                    desc = new_desc
        elif descriptor_type == 'label':
            print "Skipping label descriptor with heading {}".format(descriptor_heading)
        else:
            raise ValueError('Unrecognized descriptor type: {}'.format(descriptor_type))

        if desc is not None:
            descriptor_dict[descriptor_heading] = desc
            if save:
                desc.save()
        

    return descriptor_dict


def stringToBool(s):
    if s.lower() == 'true':
        return True
    elif s.lower() == 'false':
        return False
    else:
        raise ValueError("Tried to convert string to boolean when string was neither 'True' nor 'False' but {}".format(s))


def parse_reactions(filename, save=False, overwrite=False, val_save_cutoff=8000):
    with open(filename, 'rb') as f:
        reader = csv.DictReader(f)

        descriptor_type_dict = reader.next()

        descriptor_dict = create_descriptors(descriptor_type_dict, save=save, overwrite=overwrite)
        
        #found = []
        #failed = []

        num_vals = []
        ord_vals = []
        bool_vals = []
        cat_vals = []
        rxns_to_bulk_create = set([])
        
        for i, rxn_entry in enumerate(reader):
            ref = reimport_reactions.convert_legacy_reference(rxn_entry['XXXtitle'])
            ps = PerformedReaction.objects.filter(convertedLegacyRef=ref)
            if ps.count() == 0:
                if ref.startswith('xxx'):
                    unmunged_ref = ref + '0'
                    print '{}: Using UNMUNGED reference {} public'.format(i, unmunged_ref)
                    ps = PerformedReaction.objects.filter(convertedLegacyRef=unmunged_ref)
                else:
                    raise RuntimeError('Found {} reactions with reference {}'.format(ps.count(), ref))

            for rxn in ps:
                if rxn in rxns_to_bulk_create and save:
                    # need to save now because this reaction already has values to bulk_create
                    BoolRxnDescriptorValue.objects.bulk_create(bool_vals)
                    print "{} boolean values saved".format(len(bool_vals))
                    NumRxnDescriptorValue.objects.bulk_create(num_vals)
                    print "{} numeric values saved".format(len(num_vals))
                    OrdRxnDescriptorValue.objects.bulk_create(ord_vals)
                    print "{} ordinal values saved".format(len(ord_vals))
                    CatRxnDescriptorValue.objects.bulk_create(cat_vals)
                    print "{} categorical values saved".format(len(cat_vals))
                    bool_vals = []
                    num_vals = []
                    ord_vals = []
                    cat_vals = []
                    rxn_to_bulk_create = set([])
                rxns_to_bulk_create.add(rxn)
                for descriptor_heading, desc in descriptor_dict.items():
                    val_string = rxn_entry[descriptor_heading]
                    if val_string != '?':
                        if isinstance(desc, NumRxnDescriptor):
                            conv = float
                            val_list = num_vals
                        elif isinstance(desc, BoolRxnDescriptor):
                            conv = stringToBool
                            val_list = bool_vals
                        elif isinstance(desc, OrdRxnDescriptor):
                            conv = int
                            val_list = ord_vals
                        elif isinstance(desc, CatRxnDescriptor):
                            conv = str
                            val_list = cat_vals
                        else:
                            raise ValueError('Unrecognized descriptor type {}'.format(type(desc)))
    
                        val = desc.updateOrNewValue(rxn, conv(val_string))
                        if val is not None:
                            val_list.append(val)
                        #val.save()
                print "Created values for reaction {}, number {}".format(ref, i)

                
                if save:
                    if len(bool_vals) > val_save_cutoff:
                        BoolRxnDescriptorValue.objects.bulk_create(bool_vals)
                        print "{} boolean values saved".format(len(bool_vals))
                        bool_vals = []
                    if len(num_vals) > val_save_cutoff:
                        NumRxnDescriptorValue.objects.bulk_create(num_vals)
                        print "{} numeric values saved".format(len(num_vals))
                        num_vals = []
                    if len(ord_vals) > val_save_cutoff:
                        OrdRxnDescriptorValue.objects.bulk_create(ord_vals)
                        print "{} ordinal values saved".format(len(ord_vals))
                        ord_vals = []
                    if len(cat_vals) > val_save_cutoff:
                        CatRxnDescriptorValue.objects.bulk_create(cat_vals)
                        print "{} categorical values saved".format(len(cat_vals))
                        cat_vals = []
        if save:
            BoolRxnDescriptorValue.objects.bulk_create(bool_vals)
            print "{} boolean values saved".format(len(bool_vals))
            NumRxnDescriptorValue.objects.bulk_create(num_vals)
            print "{} numeric values saved".format(len(num_vals))
            OrdRxnDescriptorValue.objects.bulk_create(ord_vals)
            print "{} ordinal values saved".format(len(ord_vals))
            CatRxnDescriptorValue.objects.bulk_create(cat_vals)
            print "{} categorical values saved".format(len(cat_vals))

    #print "Sucessfully found {} reaction entries".format(len(found))
    #print "Failed to find {} reaction entries. Skipped.".format(len(failed))

def startswith_lookup(start):
    print list(PerformedReaction.objects.filter(reference__startswith=start))
        

    
if __name__ == '__main__':
    filename = argv[1]

    parse_reactions(filename, save=True, overwrite=False, val_save_cutoff=5000)

