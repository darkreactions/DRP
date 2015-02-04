
def apply_filters(request, data):
  """
  Reads a `request`/dictionary to parse filters and applies those filters
  on some `data` QuerySet.
  """

  from DRP.retrievalFunctions import filter_data

  filters = {}
  for key, val in request.GET.items():
    try:
      filters[key] = request.GET.getlist(key)
    except:
      filters[key] = request.GET[key]

  # If filters exist, apply them!
  if filters:
    data = filter_data(data, filters)

  return data
