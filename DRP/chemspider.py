from models.CompoundEntry import CompoundEntry
import chemspipy

def chemspider_find(search_fields):

  for i in search_fields:
    try:
      query = chemspipy.find_one(i)
      return query
    except Exception as e:
      pass
  return None


def search_chemspider(val):

  # Variable Setup
  chemspi_query = ""
  search_fields = ["CAS_ID", "compound"]

  #Accept either a CompoundEntry object or a dict with the valid fields

  if type(val)==dict:
    search_fields = [val.get(i) for i in search_fields if val.get(i)]
  elif type(val)==CompoundEntry:
    query_criteria = [val.abbrev, val.compound]
  elif type(val)==list:
    search_fields = val
  else:
    search_fields = [val]

  # `result` is a chemspider object.
  result = chemspider_find(search_fields)

  return result

