from models.CompoundEntry import CompoundEntry
import chemspipy

def chemspider_find(search_terms):
    '''Find a single chemspider instance based on search terms provided. If different search terms return a different first result this will raise a ChemspiderError'''
    #the results of this search
    results = set()
    for i in search_terms:
        result = chemspipy.find_one(i)
        resultPresent = False
        for r in results: #we shouldn't have to do this with a set, but chemspipy counts two instances of Compound with the same chemspider ID as still not equal for some reason, so there we are...
            if result.csid == r.csid:
                resultPresent = True
        if not resultPresent:
            results.add(result)

    if len(results) > 1:
        emessage = 'Conflict in results between different search terms; compounds :'
        for r in results:
            emessage += str(r)
        raise ChemspiderError(emessage)
    else:
        return results.pop()


def search_chemspider(val):
    '''Thin conformance layer for chempsider_find; can accept a list, dictionary or a CompoundEntry as an argument on which to search chemspider.'''
  # Variable Setup
    search_fields = ["CAS_ID", "compound"]
    search_terms = []

  #Accept either a CompoundEntry object or a dict with the valid fields

    if type(val)==dict:
        search_terms = [val.get(i) for i in search_fields if val.get(i)]
    elif type(val)==CompoundEntry:
        search_terms = [val.compound]
    elif type(val)==list:
        search_terms = val
    else:
        search_terms = [val]

  # `result` is a chemspider object.
    result = chemspider_find(search_terms)

    return result

class ChemspiderError(Exception):

    def __init__(s, message):
        s.message = message

    def __str__(s):
        return s.message
