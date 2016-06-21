import os
import sys
full_path = os.path.dirname(os.path.realpath(__file__)) + "/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
    sys.path = [django_path] + sys.path
    os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def get_compound_map():
    def clean_row(row):
        cleaned = []
        for elem in row:
            if elem and "\"" == elem[0] and "\"" == elem[-1]:
                elem = elem[1:-1]
            cleaned.append(elem)
        return cleaned

    import csv

    with open("raw/amines.csv") as f:
        matrix = [clean_row(line) for line in csv.reader(f)]

    headers = matrix.pop(0)
    comp_header = "CHEMICAL_NAME"
    comp_index = headers.index(comp_header)

    return {row[comp_index]: {header: row[i] for i, header in enumerate(headers)}
            for row in matrix}


def fix_wrong():
    from DRP.models import CompoundEntry
    update_refs("Diisopropyl-1,5-pentan", CompoundEntry.objects.get(id=216))


def update_refs(name, new_comp):
    from DRP.models import Recommendation

    Recommendation.objects.filter(
        reactant_fk_1__compound__icontains=name).update(reactant_fk_1=new_comp)
    Recommendation.objects.filter(
        reactant_fk_2__compound__icontains=name).update(reactant_fk_2=new_comp)
    Recommendation.objects.filter(
        reactant_fk_3__compound__icontains=name).update(reactant_fk_3=new_comp)
    Recommendation.objects.filter(
        reactant_fk_4__compound__icontains=name).update(reactant_fk_4=new_comp)
    Recommendation.objects.filter(
        reactant_fk_5__compound__icontains=name).update(reactant_fk_5=new_comp)


def update_compounds():

    from DRP.models import Recommendation
    compound_map = get_compound_map()
    compound_names = compound_map.keys()
    changed = 0

    ignore = {"piperazine", "ethylenediamine", "diethylamine", "1-methylpiperazine",
              "1,5-diaminopentane", "1,6-diaminohexane", "2-methylpiperazine",
              "1,6-hexanediamine", "4-Aminopiperidine"}

    recs = Recommendation.objects.all()
    for rec in recs:

        compounds = [rec.reactant_fk_1, rec.reactant_fk_2, rec.reactant_fk_3,
                     rec.reactant_fk_4, rec.reactant_fk_5]
        compounds = [c for c in compounds if c is not None]

        for comp in compounds:

            for name in compound_names:

                if (comp.compound.lower() in name.lower() and comp.compound != name and
                    len(comp.compound) > 5 and comp.compound not in ignore and
                        comp.compound[-3:] not in {"ine", "ane"}):

                    print "{} --> {}".format(comp.compound, name)

                    comp.compound = name
                    comp.abbrev = name

                    update_refs(name, comp)

                    comp.save()

                    changed += 1
                    break

    print "Fixed {}".format(changed)


if __name__ == "__main__":
    update_compounds()
    fix_wrong()
