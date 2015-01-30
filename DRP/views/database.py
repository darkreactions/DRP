from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
from django.shortcuts import render

@login_required
@require_http_methods(["GET"])
def database(request, page_request=None, model="database"):
  """
  The main view for interacting with Data-like elements (including Recommendations).
  """

  from DRP.models import get_lab_Data, get_recommendations
  from DRP.filters import apply_filters
  from DRP.pagifier import pagify_data, calc_total_pages, get_page_links

  # If a valid page number was given, use it.
  try:
    page = int(page_request)
  except:
    page = 1

  # Get that data that this lab group is allowed to see.
  lab = request.user.get_profile().lab_group

  if model=="database":
    data = get_lab_Data(lab)
  elif model=="recommendations":
    data = get_recommendations(lab)
  else:
    raise Exception("Invalid model specified!")

  data = apply_filters(request, data)

  # Calculate the page information.
  num_data = data.count()
  total_pages = calc_total_pages(num_data)

  # Prepare the data on a given page as (index, datum) tuples.
  data_tups = pagify_data(data, page)

  # Return a package of page information and data.
  return render(request, 'global_page.html', {
    "data": data_tups, #Includes data and data_indexes.
    "total_data_size": num_data,
    "current_page":page,
    "total_pages": total_pages,
    "page_links": get_page_links(page, total_pages),
    "query":"?"+request.GET.urlencode(),
    "template": model,
  })



