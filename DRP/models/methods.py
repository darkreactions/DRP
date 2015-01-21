
def convert_QuerySet_to_list(query, model, with_headings=True):
 #Get the appropriate headings.
 all_fields = get_model_field_names(model=model, collect_ignored=True)
 fields_to_exclude = {"lab_group", "atoms"}
 headings = [field for field in all_fields if field not in fields_to_exclude]

 if with_headings:
  query_list = [list(headings)]
 else:
  query_list = []

 for entry in query:
  sub_list = [getattr(entry, field) for field in headings]
  query_list.append(sub_list)

 return query_list


def get_model_field_names(both=False, verbose=False, model="Data",
                          unique_only=False, collect_ignored=False,
                          for_upload=False):
  clean_fields = []

  if model=="Data":
    from DRP.models.Data import Data
    all_fields = Data._meta.fields

    if collect_ignored:
      fields_to_ignore = {u"id", "creation_time_dt", "calculations"}
    else:
      fields_to_ignore = {u"id","user","lab_group", "atoms", "creation_time_dt",
                          "calculations", "calculated_temp", "calculated_time",
                          "calculated_pH", "is_valid", "public"}

  elif model=="Recommendation":
    from DRP.models.Recommendation import Recommendation
    all_fields = Recommendation._meta.fields

    if collect_ignored:
      fields_to_ignore = {u"id", "creation_time_dt"}
    else:
      fields_to_ignore = {u"id","user", "assigned_user", "lab_group", "saved",
                          "model_version", "atoms", "creation_time_dt", "nonsense",
                          "complete", "score", "date_dt", "hidden", "seed", "seeded"}
  elif model=="CompoundEntry":
    from DRP.models.CompoundEntry import CompoundEntry
    all_fields = CompoundEntry._meta.fields

    if collect_ignored:
      fields_to_ignore = {u"id", "image_url", "custom", "calculations"}
    else:
      fields_to_ignore = {u"id","lab_group", "smiles", "mw", "custom",
                          "calculations", "calculations_failed"}
  else:
    raise Exception("Unknown model specified.")


  dirty_fields = [field for field in all_fields
                                  if field.name not in fields_to_ignore]

  #Ignore any field that is in fields_to_ignore.
  for field in dirty_fields:
    #Return the non list-fields:
    if unique_only and field.name[-1].isdigit(): continue

    #Return either the verbose names or the non-verbose names.
    if both:
      clean_fields += [{"verbose":field.verbose_name, "raw":field.name}]
    elif verbose:
      clean_fields += [field.verbose_name]
    else:
      clean_fields += [field.name]

  return clean_fields

