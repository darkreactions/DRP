def smart_split(data, headers, response, reactants):
  from DRP.model_building.test_train_split import create_test_and_train_lists
  from DRP.model_building.load_data import create_reactant_keys

  response_ind = headers.index(response)

  # Create reactant-combination keys for each data entry.
  dataKeys = create_reactant_keys(data, headers=headers, reactants=reactants)

  # Partitions the data/keys into separate test/training datasets.
  test, train = create_test_and_train_lists(data, dataKeys)

  # Separate the test/training sets into predictors and responses
  A_preds = [[elem for i, elem in enumerate(row) if i!=response_ind] for row in train]
  A_resps = [row[response_ind] for row in train]
  B_preds = [[elem for i, elem in enumerate(row) if i!=response_ind] for row in test]
  B_resps = [row[response_ind] for row in test]
  splits = (A_preds, B_preds, A_resps, B_resps)

  return splits, headers
