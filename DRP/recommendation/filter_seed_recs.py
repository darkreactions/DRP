
from DRP.retrievalFunctions import *
from DRP.database_construction import *

# Removes "unusable" amines from the recList as per Alex's
#  06-10-14 request.


def filterSeedRecList(lab, recTupleList):
    # The result that will eventually be a filtered version of recTupleList.
    filteredRecs = []

    # The list of "amines that we want to be included in the ... results".
    valid_reactants = {
        "3-Aminopyrrolidine dihydrochloride",
        "N-Methylethylenediamine",
        "N,N,N',N'-Tetramethyl-1,4-butanediamine",
        "N,N,N'-Trimethyl-1,3-propanediamine",
        "4-Aminopiperidine",
        "N-Propyl-1,3-propanediamine",
        "N-Propylethylenediamine",
        "N-Isopropyl-1,3-propanediamine",
        "4-Amino-2,2,6,6-tetramethylpiperidine",
        "4-(Aminomethyl)piperidine",
        "N,N'-Diisopropylethylenediamine",
        "1-(2-Aminoethyl)piperidine",
        "N,N'-Dimethyl-1,3-propanediamine",
        "2-(Aminomethyl)-1-ethylpyrrolidine",
        "1-Butylpiperazine",
        "2-(2-Aminoethyl)-1-methylpyrrolidine",
        "N,N'-Dimethyl-1,6-hexanediamine",
        "1-(2-Aminoethyl)pyrrolidine",
        "1-Octylpiperazine",
        "4-(1-Pyrrolidinyl)piperidine",
        "2,3-Dimethylpiperazine",
        "N,N-Dimethyl-3-pyrrolidinamine",
        "3-(aminomethyl)piperidine",
        "1-(3-aminopropyl)pyrrolidine",
        "2-aminopiperidine",
        "5-Amino-1,3,3-trimethylcyclohexanemethylamine",  # TODO: add cis and trans
        "4,4'-Methylenebis(2-methylcyclohexylamine)",  # TODO: add isomers
        "N,N,N',N'-Tetramethyl-1,3-propanediamine",
        "3-(Dimethylamino)-1-propylamine",
        "1,3-Cyclohexanebis(methylamine)",  # TODO:Add isomers
        "N,N,N',N'-Tetramethyl-1,3-butanediamine",
        "2-(Aminomethyl)piperidine",
        "3-(Dimethylamino)-1-propylamine",
        "N,N,2,2-Tetramethyl-1,3-propanediamine",
        "N-Methyl-1,3-diaminopropane",
        "2-Amino-5-diethylaminopentane",
        "3-(Dibutylamino)propylamine"}

    # Get the valid Recommendation fields.
    fields = get_model_field_names(model="Recommendation")

    # The elements of recList are detailed below:
    #[(best_conf, sim(best_candidate, raw_rxn), best_candidate), ...]

    # Note that "best_candidate" is actually a reaction in the form
    #  of a list. Thus, below, we simply iterate through the reaction.

    # Also note that recList[2][0] contains the headers for the recList.

    # Actually "translate" the compounds to abbrevs.
    for recTuple in recTupleList:
        for (field, i) in zip(fields, range(len(recTuple[2][1:]))):
            if field[:8] == "reactant":
                # If there is an overlap in one of the reactants, throw it in.
                value = recTuple[2][i + 1]
                if value in valid_reactants:
                    filteredRecs.append(recTuple)
                    break

    return filteredRecs
