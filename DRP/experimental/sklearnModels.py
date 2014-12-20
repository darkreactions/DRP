
def load_data():
  import csv
  filemodel_type = "../../CSVs/12_19_14.csv"

  with open(filemodel_type) as f:
    reader = csv.reader(f)
    data = [row for row in reader]
    headers = data.pop(0)

  blacklist = {}
  cat_list = {"XXXtitle","XXXinorg1","XXXinorg2","XXXinorg3",
              "XXXorg1","XXXorg2","XXXoxlike1"}

  # Remove any columns that we don't want to use.
  i_to_remove = {i for i, header in enumerate(headers) if header in blacklist}
  headers = [h for i, h in enumerate(headers) if i not in i_to_remove]
  data = [[e for i, e in enumerate(row) if i not in i_to_remove] for row in data]

  # Categorize any non-numeric data.
  h_to_categorize = {i for i, header in enumerate(headers) if header in cat_list}
  value_map = {"yes":1, "no":0, "?":-1}
  c = 2
  for i, row in enumerate(data):
    for j, elem in enumerate(row):
      if elem in value_map:
        data[i][j] = value_map[elem]
      elif j in h_to_categorize:
        value_map[elem] = c
        data[i][j] = value_map[elem]
        c += 1

  return data, headers


def get_model(model_type):
  from sklearn.ensemble import RandomForestClassifier

  if model_type=="random forest":
    descriptors = {"n_estimators":500, "criterion":"gini"}
    model = RandomForestClassifier(**descriptors)
  else:
    raise Exception("Model model_type '{}' unknown by get_model".format(model_type))

  return model, descriptors


def split_data(data, headers, response, split=0.5):
  from sklearn.cross_validation import train_test_split

  header_index = headers.index(response)
  X = [[elem for i, elem in enumerate(row) if i!=header_index] for row in data]
  y = [row[header_index] for row in data]

  splits = train_test_split(X, y, train_size=split)
  return splits


def prepare(model, predictors, responses):
  model.fit(predictors, responses)


def test(model, predictors, responses):
  guesses = model.predict(predictors)

  possible_vals = sorted(list(set(responses))) #Assume: all values in `responses`
  len_generator = xrange(len(possible_vals))
  cm = {r1:{r2:0 for r2 in possible_vals} for r1 in possible_vals}

  for guess, response in zip(guesses, responses):
    cm[response][guess] += 1

  return cm


def analyze(model):
  model.graph_confusion_table()
  model.pretty_stats()


def make_sklearn_ModelStats(sklearn_model, cm, response, title="", description=""):
  from DRP.models import ModelStats
  model = ModelStats()
  model.autofill(title=title, description=description)

  model.set_confusion_table(cm)
  model.set_correct_vals(["3.0","4.0"])

  model.save_model_file(sklearn_model, use_sklearn=True)

  return model


def main():
  model_type = "random forest"
  response = "outcome"
  train_percentage = 0.5

  # Organize and split the data.
  data, headers = load_data()
  sklearn_model, descriptors = get_model(model_type)
  A_preds, B_preds, A_resps, B_resps = split_data(data, headers, response,
                                                  split=train_percentage)

  # Build and test the sklearn model.
  prepare(sklearn_model, A_preds, A_resps)
  cm = test(sklearn_model, B_preds, B_resps)

  # Construct and store the model database entry.
  descriptors["response"] = response
  descriptors["model_type"]=model_type
  descriptors["train percentage"]=train_percentage
  title = "'{}' on '{}'".format(model_type, response)
  description = ", ".join(["{}:{}".format(key,val) for key,val in descriptors.items()])

  model = make_sklearn_ModelStats(sklearn_model, cm, A_resps,
                                  title=title, description=description)

  analyze(model)

