
import sys, os
django_dir = os.path.dirname(os.path.realpath(__file__)).split("DRP")[0]
django_path = "{}/DRP".format(django_dir)
if django_path not in sys.path:
  sys.path.append("{}/DRP".format(django_dir))

os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


def default_splitter(data, headers=None):
  from DRP.model_building.test_train_split import create_test_and_train_lists
  from DRP.model_building.load_data import create_reactant_keys

  # Create reactant-combination keys for each data entry.
  dataKeys = create_reactant_keys(data, headers=headers)

  # Partitions the data/keys into separate test/training datasets.
  test, train = create_test_and_train_lists(data, dataKeys)

  return {"test":test, "train":train, "all":data}


def naive_split(data, headers=None):
  import random

  test_fraction = 0.3
  split_index = int(len(data)*test_fraction)

  random.shuffle(data)

  return{"all":data, "test":data[:split_index], "train":data[split_index:]}





