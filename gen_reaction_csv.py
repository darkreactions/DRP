import csv

from DRP.models import PerformedReaction

rxns = PerformedReaction.objects.all().rows(expanded=True)

#print (type(rxns))

key_set = set()
for i in rxns:
    for key in i:
        key_set.add(key)

sorted_keys = sorted(key_set)

# Write to csv
with open('Database_Reactions.csv', 'w') as f:
    writer = csv.DictWriter(f, fieldnames = sorted_keys)
    rxns = PerformedReaction.objects.all().rows(expanded=True)
    writer.writeheader()
    for rxn in rxns:
        writer.writerow(rxn)

    
