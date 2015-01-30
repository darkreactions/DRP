from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
from django.shortcuts import render

@login_required
@require_http_methods(["GET"])
def database(request, page_request=None):
  from DRP.models import get_lab_Data
  from DRP.retrievalFunctions import filter_data
  from DRP.pagifier import pagify_data, calc_total_pages, get_page_links

  # If a valid page number was given, use it.
  try:
    page = int(page_request)
  except:
    page = 1

  # Get that data that this lab group is allowed to see.
  lab = request.user.get_profile().lab_group
  data = get_lab_Data(lab)

  filters = {}
  for key, val in request.GET.items():
    try:
      filters[key] = request.GET.getlist(key)
    except:
      filters[key] = request.GET[key]

  # If filters exist, apply them!
  if filters:
    data = filter_data(data, filters)

  # Calculate the page information.
  num_data = data.count()
  total_pages = calc_total_pages(num_data)

  # Prepare the data on a given page.
  data = pagify_data(data, page)

  from DRP.data_config import CONFIG
  start_index = (page-1) * CONFIG.data_per_page
  data_tups = [(elem, start_index+i+1) for i, elem in enumerate(data)]

  # Return a package of page information and data.
  return render(request, 'global_page.html', {
    "data": data_tups, #Includes data and data_indexes.
    "total_data_size": num_data,
    "current_page":page,
    "total_pages": total_pages,
    "page_links": get_page_links(page, total_pages),
    "query":"?"+request.GET.urlencode(),
    "template": "database",
  })



