
import sys, os
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def default_splitter(data, headers=None):
  from DRP.model_building.test_train_split import create_test_and_train_lists
  from DRP.model_building.load_data import create_reactant_keys
  import random

  random.shuffle(data)

  # Create reactant-combination keys for each data entry.
  dataKeys = create_reactant_keys(data, headers=headers)

  # Partitions the data/keys into separate test/training datasets.
  test, train = create_test_and_train_lists(data, dataKeys)

  return {"test":test, "train":train, "all":data}


def add_nonsense_to_test(data, headers=None):
  from DRP.models import Recommendation
  from DRP.preprocessors import default_preprocessor

  splits = default_splitter(data, headers=headers)

  print "Calculating test supplement..."
  nonsense_recs = Recommendation.objects.filter(nonsense=True)
  nonsense = default_preprocessor([headers]+list(nonsense_recs))
  nonsense = [row[:-1] + [0] for row in nonsense]
  print "Supplement size: {}".format(len(nonsense))


  splits["test"].extend(nonsense[1:])
  splits["all"].extend(nonsense[1:])

  return splits




def naive_splitter(data, headers=None):
  import random

  test_fraction = 0.3
  split_index = int(len(data)*test_fraction)

  random.shuffle(data)

  return {"all":data, "test":data[:split_index], "train":data[split_index:]}

def strict_category_splitter(data, headers=None):

  # If the last field in a datum is `True` let it be in the `test` set.

  test, train = [], []
  for datum in data:
    if datum.pop(-1)==True:
      test.append(datum)

    else:
      train.append(datum)

  return {"test":test, "train":train, "all":data}


def category_splitter(data, headers=None):

  # Remove any entries from the "test" data that are not in the category.

  splits = default_splitter(data, headers=headers)

  non_category_test = []
  category_test = []

  for datum in splits["test"]:
    if datum[-1]==True:
      category_test.append(datum)
    else:
      non_category_test.append(datum)


  remove_last = lambda row: row[:-1]

  splits = {
    "all": map(remove_last, splits["all"]),
    "train": map(remove_last, splits["train"]+non_category_test),
    "test": map(remove_last, category_test),
  }

  return splits

