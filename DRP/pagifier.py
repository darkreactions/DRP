def get_page_links(current, total):
  if not total: return []

  from data_config import CONFIG

  #Variable Setup:
  radius = CONFIG.current_page_radius

  #Always display at least page one.
  page_set = {1}
  if total==1:
    return list(page_set)
  else:
    for i in xrange(current-radius, current+radius+1):
      if (1 < i < total): page_set.add(i)

    page_set.add(total) #Always show the last page.

    #Convert page_links to an ordered list.
    raw_page_links = list(page_set)
    raw_page_links.sort()

    page_links = []
    for i in raw_page_links:
      if page_links and i-page_links[-1]>1:
        #Add an ellipsis between any pages that have a gap between them.
        page_links.append("...")
        page_links.append(i)
      else:
        page_links.append(i)
    return page_links



def calc_total_pages(data_size):
  from DRP.data_config import CONFIG
  total_pages = 1 + int((data_size-1)/CONFIG.data_per_page)
  return total_pages if total_pages > 1 else 1



def pagify_data(data, page):
    from DRP.data_config import CONFIG

    #Variable Setup:
    data_per_page = CONFIG.data_per_page

    total_pages = calc_total_pages(data.count())

    #Check that the page request is valid.
    if not (0<page<=total_pages):
        return []

    #Return the data that would be contained on the requested page.
    page_data = data[(page-1)*data_per_page:page*data_per_page]
    return page_data



def repackage_page_session(session):
  from DRP.data_config import CONFIG
  data = session["page_data"]
  total_pages = session["total_pages"]
  page_links = session["page_links"]
  current_page = session["current_page"]
  total_data_size = session["total_data_size"]

  #Show the overall index of each datum.
  start_index = (current_page-1)*CONFIG.data_per_page + 1
  end_index = (current_page)*CONFIG.data_per_page + 1

  #Prepare packages.
  data_package = zip(data, range(start_index, end_index))
  page_info = {
    "current_page":current_page,
    "total_pages":total_pages,
    "data_per_page":CONFIG.data_per_page,
    "page_links":page_links,
    }
  return data_package, page_info, total_data_size

