def default_postprocessor(splits, headers):
  """
  The default post-processor to be used in model-generation
  after the data-splitting process. Removes unnecessary headers
  and sets values to those expected by the model.
  """

  # Remove unnecessary fields from the data.
  blacklist_fields = [
  "XXXtitle", "XXXinorg1", "XXXinorg2", "XXXinorg3",
  "XXXorg1", "XXXorg2", "XXXoxlike1", "XXXinorg1mass",
  "XXXinorg1moles", "XXXinorg2mass", "XXXinorg2moles",
  "XXXinorg3mass", "XXXinorg3moles", "XXXorg1mass",
  "XXXorg1moles", "XXXorg2mass", "XXXorg2moles",
  "XXXoxlike1mass", "XXXoxlike1moles", "purity",
  ]

  blacklist = {headers.index(field)
                 for field in blacklist_fields if field in headers}

  headers = [field for i, field in enumerate(headers) if i not in blacklist]
  splits = { key:[ [e for i,e in enumerate(row) if i not in blacklist]
                        for row in data]
                          for key, data in splits.items()}


  # Convert any stringy-bools in the data to integers.
  #TODO: Erase
  """
  bool_map = {
  "yes":1,
  "no":0,
  "?":-1
  }

  for key, data in splits.items():
    for i, row in enumerate(data):
      for j, elem in enumerate(row):
        if elem in bool_map:
          splits[key][i][j] = bool_map[elem]
  """

  return (splits, headers)


def reduce_features(splits, headers):
  from DRP.model_building.rxn_calculator import josh_serota_feature_selected

  headers_to_keep = josh_serota_feature_selected
  indexes = {headers.index(header) for header in headers_to_keep}

  # Remove the dimensions that are not selected.
  headers = [elem for i, elem in enumerate(headers) if i in indexes]

  splits = { key:[ [e for i,e in enumerate(row) if i in indexes]
                        for row in data]
                          for key, data in splits.items()}

  return default_postprocessor(splits, headers)



