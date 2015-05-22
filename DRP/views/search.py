from django.shortcuts import render


def search(request, model="Data", params={}):
  from DRP.models import get_model_field_names

  search_fields = get_model_field_names(both=True, unique_only=True, model=model)

  if model=="Data":
    #Collect the fields that will be displayed in the Search "Fields" tab.
    search_fields = get_model_field_names(both=True, unique_only=True)
    search_fields = [
    {"raw":"reactant", "verbose":"Reactant"},
    {"raw":"quantity", "verbose":"Quantity"},
    {"raw":"unit", "verbose":"Unit"},
    {"raw":"is_valid", "verbose":"Is Valid"},
    {"raw":"user", "verbose":"User"},
    {"raw":"public", "verbose":"Public"}] + search_fields
  elif model=="Recommendation":
    #Collect the fields that will be displayed in the Search "Fields" tab.
    search_fields = [
    {"raw":"reactant", "verbose":"Reactant"},
    {"raw":"quantity", "verbose":"Quantity"},
    {"raw":"unit", "verbose":"Unit"},
    {"raw":"assigned_user", "verbose":"Assigned User"}] + search_fields

    if "seeded" in params and params["seeded"]:
      model = "SeedRecommendation"
      search_fields += [{"raw":"seed", "verbose":"Seed is..."}]


  search_sub_fields = [
    {"raw":"abbrev", "verbose":"Abbrev"},
    {"raw":"compound", "verbose":"Compound"},
    {"raw":"compound_type", "verbose":"Type"}
  ]

  return render(request, 'search_global.html', {
   "search_fields": search_fields,
   "search_sub_fields": search_sub_fields,
   "model": model,
  })

