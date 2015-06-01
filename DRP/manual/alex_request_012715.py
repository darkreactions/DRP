
import os, sys
full_path = os.path.dirname(os.path.realpath(__file__))+"/"
django_path = full_path[:full_path.rfind("/DRP/")]
if django_path not in sys.path:
    sys.path = [django_path] + sys.path
    os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def get_page(datum):
    import re
    ref = datum.ref
    cleaned = re.sub(r"[^\d\.]", "", ref)
    parts = cleaned.split(".")

    if len(parts[1]) == 1:
      parts[1] = "0" + parts[1]

    return float("{0[0]}.{0[1]}".format(parts))


def strip(datum, headers):
    from DRP.models import CompoundEntry

    row = []

    for header in headers:
        if "." in header:
            model, field = header.split(".")
            if "reactant" in model:
                abbrev = getattr(datum, model)
                foreign_key = CompoundEntry.objects.get(abbrev=abbrev)
            else:
                foreign_key = getattr(datum, model)
            val = getattr(foreign_key, field)
        else:
            val = getattr(datum, header)

        row.append(val)

    return row


def write_csv(name, headers, matrix):
    import csv

    with open(name, "w") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in matrix:
            writer.writerow(row)



def main():

    from DRP.models import Data

    groups = {i:[] for i in xrange(1,4)}

    for elem in Data.objects.all():

        if "ke" in elem.ref:
            page = get_page(elem)

            if 17.01 <= page <= 17.12:
                groups[1].append(elem)
            elif 22.01 <= page <= 22.27:
                groups[1].append(elem)

            elif 28.01 <= page <= 28.39:
                groups[2].append(elem)
            elif 61.01 <= page <= 61.13:
                groups[2].append(elem)
            elif 70.01 <= page <= 70.25:
                groups[2].append(elem)

            elif 56.01 <= page <= 56.25:
                groups[3].append(elem)
            elif 64.01 <= page <= 64.13:
                groups[3].append(elem)

    headers = ["ref", "reactant_3.abbrev", "reactant_3.compound",
               "outcome", "reactant_3.smiles"]

    for group, elems in groups.items():

        new_elems = [strip(elem, headers) for elem in elems]

        write_csv("group_{}.csv".format(group), headers, new_elems)


    print "DONE!"

if __name__=="__main__":
    main()
