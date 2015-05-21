
def apply_filters(request, queryset, model="Data"):
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
    if model=="Data":
      queryset = filter_data(queryset, filters)
    elif model=="ModelStats":
      queryset = filter_models(queryset, filters)

  return queryset


def filter_models(models, filters):
  from django.db.models import Q
  import operator

  queries = []
  for key, vals in filters.items():
    Qs = [Q(**{key:val}) for val in vals]
    query = reduce(operator.or_, Qs)
    queries.append(query)

  overall_query = reduce(operator.and_, queries)
  models = models.filter(overall_query)
  return models


