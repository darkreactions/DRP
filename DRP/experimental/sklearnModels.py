
def load_data():
  import csv
  filemodel_type = "../../CSVs/12_19_14.csv"

  with open(filemodel_type) as f:
    reader = csv.reader(f)
    data = [row for row in reader]
    headers = data.pop(0)

  # Remove any columns that we don't want to use.
  blacklist = {}
  i_to_remove = {i for i, header in enumerate(headers) if header in blacklist}
  headers = [h for i, h in enumerate(headers) if i not in i_to_remove]
  data = [[e for i, e in enumerate(row) if i not in i_to_remove] for row in data]

  # Categorize any non-numeric data.
  c = 2
  value_map = {"yes":1, "no":0}
  for i, row in enumerate(data):
    for j, elem in enumerate(row):
      if elem in value_map:
        data[i][j] = value_map[elem]

  return data, headers


def get_model(model_type):
  if model_type=="random forest":
    from sklearn.ensemble import RandomForestClassifier
    model = RandomForestClassifier(n_estimators=100, criterion="gini")
  else:
    raise Exception("Model model_type '{}' unknown by get_model".format(model_type))

  return model


def split_data(data, headers, response):
  from sklearn.cross_validation import train_test_split
  test_size = 0.7

  header_index = headers.index(response)
  X = [[elem for i, elem in enumerate(row) if i!=header_index] for row in data]
  y = [row[header_index] for row in data]

  splits = train_test_split(X, y, test_size=test_size)
  return splits


def prepare(model, predictors, responses):
  model.fit(predictors, responses)


def test(model, predictors, responses):
  guesses = model.predict(predictors)

  possible_vals = sorted(list(set(responses))) #Assume: all values in `responses`
  len_generator = xrange(len(possible_vals))
  cm = [[0 for j in len_generator] for i in len_generator]

  for guess, response in zip(guesses, responses):
    g_index = possible_vals.index(guess)
    r_index = possible_vals.index(response)
    cm[r_index][g_index] += 1

  return cm


def main():
  model_type = "random forest"
  response = "outcome"

  data, headers = load_data()
  model = get_model(model_type)
  X_train, y_train, X_test, y_test = split_data(data, headers, response)
  prepare(model, X_train, y_train)
  cm = test(model, X_test, y_test)

  print cm

