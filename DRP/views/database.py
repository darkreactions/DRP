from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
from django.shortcuts import render, redirect

@login_required
@require_http_methods(["GET"])
def database(request, page_request=None, model="database"):
  """
  The main view for interacting with Data-like elements (including Recommendations).
  """

  from DRP.models import get_lab_Data, get_recommendations, get_users
  from DRP.filters import apply_filters
  from DRP.pagifier import pagify_data, calc_total_pages, get_page_links

  # If a valid page number was given, use it.
  try:
    page = int(page_request)
  except:
    page = 1

  formatting_info = {}
  current_query = "?"+request.GET.urlencode(),

  # Get that data that this lab group is allowed to see.
  lab = request.user.get_profile().lab_group

  if model=="database":
    data = get_lab_Data(lab)
    formatting_info["show_logistics"] = True
    formatting_info["copyable"] = True
    formatting_info["allow_uploads"] = False

  elif model=="recommendations":
    data = get_recommendations(lab)
    formatting_info["recommendations"] = True

  elif model=="saved":
    data = get_recommendations(lab).filter(saved=True)
    formatting_info["transferable"] = True
    formatting_info["recommendations"] = True
    formatting_info["users"] = get_users(lab)

  else:
    raise Exception("Invalid model specified!")

  data = apply_filters(request, data)

  # Calculate the page information.
  num_data = data.count()
  total_pages = calc_total_pages(num_data)

  # Prepare the data on a given page as (index, datum) tuples.
  data_tups = pagify_data(data, page)

  if not data_tups and page>1:
    new_url = "/database/"
    return redirect(new_url)


  # Return a package of page information and data.
  return render(request, 'global_page.html', {
    "data": data_tups, #Includes data and data_indexes.
    "total_data_size": num_data,
    "current_page":page,
    "total_pages": total_pages,
    "page_links": get_page_links(page, total_pages),
    "query":current_query,
    "model":model,
    "template": "database",
    "formatting":formatting_info,
  })



