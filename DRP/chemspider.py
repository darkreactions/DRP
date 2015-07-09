from models.CompoundEntry import CompoundEntry
import chemspipy

def chemspider_find(search_terms):

  for i in search_terms:
      query = chemspipy.find_one(i)
      return query


def search_chemspider(val):

  # Variable Setup
  search_fields = ["CAS_ID", "compound"]
  search_terms = []

  #Accept either a CompoundEntry object or a dict with the valid fields

  if type(val)==dict:
    search_terms = [val.get(i) for i in search_fields if val.get(i)]
  elif type(val)==CompoundEntry:
    query_terms = [val.compound]
  elif type(val)==list:
    search_terms = val
  else:
    search_terms = [val]

  # `result` is a chemspider object.
  result = chemspider_find(search_terms)

  return result

