def get_page_link_format(current, total):
  if not total: return []

  from data_config import CONFIG

  #Variable Setup:
  radius = CONFIG.current_page_radius
  data_per_page = CONFIG.data_per_page

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

