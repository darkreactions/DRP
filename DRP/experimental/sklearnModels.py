def clean_data(data, headers):
  cloned_data = [row[:] for row in data]

  outcome_index = headers.index("outcome")

  for i, row in enumerate(data):
    cloned_data[i][outcome_index] = row[outcome_index] if row[outcome_index] else 1

  return cloned_data


def load_data(csv_to_use=""):
  import csv
  from DRP.model_building.load_data import load
  from DRP.retrievalFunctions import get_valid_data,expand_data,get_expanded_headers

  # Variable Setup
  categorize_nonnumeric = True

  if csv_to_use:
    with open(csv_to_use) as f:
      reader = csv.reader(f)
      data = [row for row in reader]
      headers = data.pop(0)
  else:
    db_data = get_valid_data()
    data = expand_data(db_data)
    headers = get_expanded_headers()

  blacklist = {}
  cat_list = {"XXXtitle","XXXinorg1","XXXinorg2","XXXinorg3",
              "XXXorg1","XXXorg2","XXXoxlike1"}

  # Remove any columns that we don't want to use.
  i_to_remove = {i for i, header in enumerate(headers) if header in blacklist}
  headers = [h for i, h in enumerate(headers) if i not in i_to_remove]
  data = [[e for i, e in enumerate(row) if i not in i_to_remove] for row in data]

  if categorize_nonnumeric:
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

  data = clean_data(data, headers)

  return data, headers

def split_data_paul(data, headers, response, split=0.5, PCs_to_use=0):
  # Calls Paul's splitter; need to implement split size and PCA.
  from DRP.experimental.splitter import smart_split
  reactants = ["XXXinorg1","XXXinorg2","XXXinorg3", "XXXorg1","XXXorg2","XXXoxlike1"]
  return smart_split(data, headers, response, reactants)


def split_data(data, headers, response, split=0.5, PCs_to_use=0):
  from sklearn.cross_validation import train_test_split
  from sklearn.decomposition import PCA

  header_index = headers.index(response)
  X = [[elem for i, elem in enumerate(row) if i!=header_index] for row in data]
  y = [row[header_index] for row in data]

  if PCs_to_use>0:
    pca = PCA(n_components=PCs_to_use)
    X = pca.fit_transform(X)
    updated_headers = ["PC{}" for i in xrange(PCs_to_use)]
  else:
    updated_headers = headers

  splits = train_test_split(X, y, train_size=split)
  return splits, updated_headers


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


def get_model(model_type):

  if model_type=="random forest":
    from sklearn.ensemble import RandomForestClassifier as model
    descriptors = {"n_estimators":500, "criterion":"gini", "n_jobs":-1}

  elif model_type=="linear regression": #TODO: "Cannot perform reduce with flexible type"
    from sklearn.linear_model import LinearRegression as model
    descriptors = dict()

  elif model_type=="svc":
    from sklearn.svm import SVC as model
    descriptors = {"C":1, "kernel":"linear"}

  elif model_type=="knn":
    from sklearn.neighbors import KNeighborsClassifier as model
    descriptors = {"p":3, "n_neighbors":1, "weights":"distance"}

  else:
    raise Exception("Model model_type '{}' unknown by get_model".format(model_type))

  return model(**descriptors), descriptors


def gen_sklearn_model(model_type, splits, headers, details={}, title=""):
  import time

  # Variable Setup
  A_preds, B_preds, A_resps, B_resps = splits

  # Build the sklearn model.
  start_time = time.time()
  sklearn_model, descriptors = get_model(model_type)
  prepare(sklearn_model, A_preds, A_resps)
  descriptors["generation time"]=time.time()-start_time

  start_time = time.time()
  cm = test(sklearn_model, B_preds, B_resps)
  descriptors["test time"]=time.time()-start_time

  # Construct and store the model database entry.
  descriptors.update(details)
  description = ", ".join(["{}:{}".format(key,val) for key,val in descriptors.items()])

  if not title:
    title = "{}_sklearn".format(str(int(time.time())))

  model = make_sklearn_ModelStats(sklearn_model, cm, A_resps,
                                  title=title, description=description)
  return model, cm

def main():
  response = "outcome"
  train_percentage = 0.5
  model_type = "random forest"
  PCs_to_use = 0

  details = {
    "model_type":model_type,
    "response":response,
    "PCs_to_use":PCs_to_use,
    "train percentage":train_percentage
  }

  title = "'{}' on '{}'".format(model_type, response)

  # Prepare the data to use in the model.
  data, headers = load_data()
  splits, headers = split_data(data, headers, response, split=train_percentage,
                                                        PCs_to_use=PCs_to_use)

  # Generate the model and store it in the database.
  model, cm = gen_sklearn_model(model_type, splits, headers, details=details, title=title)

  analyze(model)
