import django
django.setup()
from itertools import chain
from DRP.models import NumRxnDescriptor, NumRxnDescriptorValue, PerformedReaction


def create_descs():
    desc_types = [NumRxnDescriptor]

    for desc_type in desc_types:
        descs_to_save = []
        fields = [field.name for field in desc_type._meta.get_fields() if (
            not field.is_relation and not field.auto_created)]
        for d_fields in desc_type.objects.filter(heading__contains='_pH1_').exclude(heading__startswith='_').values(*fields):
            d_fields['heading'] = d_fields[
                'heading'].replace('pH1', 'pHreaction')
            d_fields['name'] = d_fields['name'].replace(
                'at pH 1', 'at reaction pH')
            d = NumRxnDescriptor(**d_fields)
            d.save()


def create_vals(rxns):
    vals_to_create = []
    count = 1
    for rxn in rxns:
        print "{} ({}/{})".format(rxn.reference, count, rxns.count())
        reaction_pH = NumRxnDescriptorValue.objects.get(
            reaction=rxn, descriptor__heading='reaction_pH').value
        if reaction_pH is not None:
            rounded_pH = int(round(reaction_pH))
            if rounded_pH == 0:
                rounded_pH = 1
            if rounded_pH not in range(1, 15):
                raise ValueError(
                    "Found pH {}. Doesn't round to [0,...,14]".format(reaction_pH))

            for d in NumRxnDescriptor.objects.filter(heading__contains='_pHreaction_'):
                headings = [d.heading.replace(
                    '_pHreaction_', '_pH{}_'.format(n)) for n in range(1, 15)]
                pH_descriptor_values = NumRxnDescriptorValue.objects.filter(
                    descriptor__heading__in=headings, reaction=rxn)
                if pH_descriptor_values.count() == 14:
                    reaction_pH_descriptor_heading = d.heading.replace(
                        '_pHreaction_', '_pH{}_'.format(rounded_pH))
                    reaction_pH_descriptor = NumRxnDescriptor.objects.get(
                        heading=reaction_pH_descriptor_heading)
                    reaction_pH_descriptor_value = pH_descriptor_values.values_list(
                        'value', flat=True).get(descriptor=reaction_pH_descriptor)
                    n = NumRxnDescriptorValue(
                        descriptor=d, reaction=rxn, value=reaction_pH_descriptor_value)
                    vals_to_create.append(n)
                    # n.save()
                    # print "val saved"
                elif pH_descriptor_values.count() != 0:
                    raise RuntimeError("Reaction {} has {} values defined for descriptor starting with {} and ending with {}. Should be 0 or 14".format(
                        rxn.reference, pH_descriptor_values.count(), prefix, suffix))

        if len(vals_to_create) > 5000:
            print "Creating values"
            NumRxnDescriptorValue.objects.bulk_create(vals_to_create)
            vals_to_create = []

        count += 1

    print "Creating values"
    NumRxnDescriptorValue.objects.bulk_create(vals_to_create)

if __name__ == '__main__':
    # create_descs()
    rxns = PerformedReaction.objects.filter(valid=True).exclude(compounds=None)

    create_vals(rxns)

    print NumRxnDescriptorValue.objects.filter(descriptor__heading__contains='_pHreaction_').count()
    # print
    # NumRxnDescriptor.objects.filter(heading__contains='_pHreaction_').count()
