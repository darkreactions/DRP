import sys
import os

django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
    sys.path.append("{}/DRP".format(django_dir))
os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'

import DRP.models
from DRP.recommendation.metrics import get_default_metric
from DRP.model_building import load_cg
from DRP.logPrinting import print_error

# Variable Setup.
CG = load_cg.get_cg()

reactant_fields = {
    "reactant_fk_1": ("quantity_1", "unit_1"),
    "reactant_fk_2": ("quantity_2", "unit_2"),
    "reactant_fk_3": ("quantity_3", "unit_3"),
    "reactant_fk_4": ("quantity_4", "unit_4"),
    "reactant_fk_5": ("quantity_5", "unit_5"),
}


pH_range = [1, 3, 5]
number_of_amine_mole_steps = 10
steps_per_amine = len(pH_range) * number_of_amine_mole_steps


def get_reactants_indices(reaction):
    indices = dict()
    for r_f, params in reactant_fields.items():
        compound = getattr(reaction, r_f)
        mass = getattr(reaction, params[0])
        units = getattr(reaction, params[1])

        # Validate this reactant...
        if compound.compound not in CG:
            continue
        if units.lower() != "g":
            raise Exception("Units not 'g'!")
        if not mass:
            raise Exception("Mass was either 0 or non-numeric!")

        # Figure out what reactants are where in the `reaction`.
        print compound, compound.compound_type
        if compound.compound_type == "Org":
            indices["org"] = r_f
        elif compound.compound_type == "Inorg":
            if "metal_1" in indices:
                indices["metal_2"] = r_f
            else:
                indices["metal_1"] = r_f

        elif compound.compound_type == "Water":
            indices["water"] = r_f

    if ("metal_1" not in indices or "metal_2" not in indices or
            "org" not in indices or "water" not in indices):
        raise Exception("Missing a crucial reactant type: " + str(indices))

    return indices


def get_reactant_moles(reaction, reactant_index):
    molar_mass = CG[getattr(reaction, reactant_index).compound]["mw"]
    unit_field = reactant_index.replace("reactant_fk", "quantity")
    mass = getattr(reaction, unit_field)
    return molar_mass / float(mass)


def get_amine_range(moles):
    return [moles * i / 1000.0 for i in xrange(0, 10)]


def row_generator(reaction, indices, amine_moles, amine_list):

    metal_1 = getattr(reaction, indices["metal_1"])
    metal_2 = getattr(reaction, indices["metal_2"])
    amine = getattr(reaction, indices["org"])
    water = getattr(reaction, indices["water"])

    metal_1_mass = getattr(reaction, reactant_fields[indices["metal_1"]][0])
    metal_2_mass = getattr(reaction, reactant_fields[indices["metal_2"]][0])
    water_mass = getattr(reaction, reactant_fields[indices["water"]][0])

    amine_range = get_amine_range(amine_moles)

    pH_range = [1, 3, 5]

    for amine in amine_list:
        if amine not in CG:
            print "Not in CG: {0}".format(amine)
            continue
        for pH in pH_range:
            for mass in amine_range:
                yield ["--", metal_1, metal_1_mass, "g", metal_2, metal_2_mass,
                       "g", amine, mass, "g", water, water_mass, "g", "", "", "",
                       getattr(reaction, "temp"), getattr(reaction, "time"), pH, "yes", "no", 4, 2, ""]

# TODO: Look here for where to change the amine molar mass, Casey.


def get_candidates(results, idx, raw_rows):
    print_error("{} results to check...".format(len(results)))
    candidates = []
    for i in range(steps_per_amine):
        if i + idx in results:
            candidates.append((results[i + idx], raw_rows[i + idx]))
    if len(candidates) == 0:
        return None, 0.0

    num_candidates = len(candidates)
    best_conf = max(candidates, key=lambda x: x[0])[0]
    candidates = filter(lambda x: x[0] == best_conf, candidates)
    candidates = [c[1] for c in candidates]

    avg_mass = 0.0
    avg_ph = 0.0
    for r in candidates:
        avg_mass += float(r[8])
        avg_ph += float(r[18])

    avg_mass = avg_mass / float(len(candidates))
    avg_ph = avg_ph / float(len(candidates))

    smallest_dist = float("inf")
    best_row = None

    for row in candidates:
        dist = (avg_mass - float(row[8]))**2 + (avg_ph - float(r[18]))**2
        if dist < smallest_dist:
            smallest_dist = dist
            best_row = row

    assert(best_row is not None)

    score = len(candidates) / float(num_candidates)

    return best_row, score


def generate_grid(reaction, amine_list, debug=True):
    from DRP.retrievalFunctions import get_latest_ModelStats

    indices = get_reactants_indices(reaction)
    amine_moles = get_reactant_moles(reaction, indices["org"])

    row_gen = row_generator(reaction, indices, amine_moles, amine_list)
    predictors = [row for row in row_gen]

    model = get_latest_ModelStats()

    if model is None:
        raise Exception("Latest model could not be found!")

    guesses = model.predict(predictors, table="tmp")

    correct = model.load_correct_vals()

    successes = [rxn for i, rxn in enumerate(
        predictors) if guesses[i] in correct]
    # TODO: NEED to grab "confidences" from the model.

    amines_results = []
    raw_rxn = DRP.models.convert_Data_to_list(reaction)
    for i, amine in enumerate(amine_list):
        best_candidate, best_conf = get_candidates(successes,
                                                   i * steps_per_amine,
                                                   predictors)
        if best_candidate:
            amine_tuple = (best_conf,
                           get_default_metric(best_candidate, raw_rxn),
                           best_candidate)
            amines_results.append(amine_tuple)

    if debug:
        print "found {} possible recommendations.".format(len(amines_results))

    amines_results.sort(key=lambda x: x[0] * x[1], reverse=True)

    return amines_results


def constructRecsFromSeed(seed_pid):
    from DRP.models import Lab_Group, CompoundEntry, Data
    # Load the default amines.
    try:
        amine_lab = Lab_Group.objects.filter(lab_title="default_amines")[0]
        amines_raw = CompoundEntry.objects.filter(lab_group=amine_lab)
        amines_names = [c.compound for c in amines_raw]
    except Exception as e:
        raise Exception("Could not load default_amines...\n{}".format(e))

    # Get the datum and lab group
    try:
        rxn = Data.objects.get(id=seed_pid)
    except Exception as e:
        raise Exception(
            "Could not use Datum \"{}\" as seed...\n{}".format(seed_pid, e))

    # Actually create the recommendations from the supplied amines and Datum.
    try:
        recommendation_list = generate_grid(rxn, amines_names, debug=False)
    except Exception as e:
        raise Exception(
            "Could not generate_grid for Datum: {}\n{}".format(seed_pid, e))

    if not recommendation_list:
        raise Exception(
            "No recommendations generated from Datum: {}".format(seed_pid))

    # Only return a number of the possible recommendations.
    return recommendation_list


# Generate seeds from a Datum's id.
if __name__ == "__main__":
    try:
        print constructRecsFromSeed(sys.argv[1])
        print "Seed Recommendations Built Successfully!"
    except Exception as e:
        print e
        print "\nUSAGE: python this_script.py <datum.id>"
