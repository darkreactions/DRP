import random

TEST_PERCENT = .30
MAX_PARTITION_SIZE = 35 # MAGIC!



# JUST KIDDING.
# I chose to only include relatively small partitions in the test set.
# Choosing 35 makes more than 95% of partitions eligible, and more than 80%
# of the data eligible.

def create_key(line):
	key = [line[0], line[3], line[6], line[9], line[12]]
	key = [r for r in key if r.lower() != 'water' and r != ""]
	key.sort()
	return tuple(key)

def create_key_in_test_map(data_list):

	number_of_reactions = 0
	key_count_map = dict()

	for line in data_list:
		key = create_key(line)
		if key not in key_count_map:
			key_count_map[key] = 0
		key_count_map[key] += 1
		number_of_reactions += 1

	key_list = key_count_map.keys()
	### Shuffle the list. We can now pick items from the top of the list
	random.shuffle(key_list)
	key_in_test = dict()
	test_total = 0
	max_test = int(TEST_PERCENT*float(number_of_reactions))

	i = 0
	for key in key_list:
		if test_total < max_test and key_count_map[key] <= MAX_PARTITION_SIZE:
			key_in_test[key] = True
			test_total += key_count_map[key]
		else:
			key_in_test[key] = False
	return key_in_test

def build_key_in_test_map(keys):
  key_counts = {k: keys.count(k) for k in set(keys)}
  key_list = key_counts.keys()
  key_in_test = dict()
  test_total = 0
  max_test = int(TEST_PERCENT*len(keys))
  i = 0
  for key in key_list:
    if test_total < max_test and key_counts[key] <= MAX_PARTITION_SIZE:
      key_in_test[key] = True
      test_total += key_counts[key]
    else:
      key_in_test[key] = False
  return key_in_test


# Returns two separate lists of data entries (a test and a training list)
#   given data entries and the corresponding keys (in the same order).
def create_test_and_train_lists(data, keys):
  test, train = [], []

  # Choose which key should go where (ie: whether a datum will
  #   be thrown to the test or training lists).
  key_in_test_map = build_key_in_test_map(keys)

  for i in xrange(len(keys)):
    datum = data[i]
    key = keys[i]

    if key_in_test_map[key]:
      test.append(datum)
    else:
      train.append(datum)

  return test, train



def is_match(line, m):
	ts = [line[1], line[4], line[7], line[10], line[13]]
	for t in ts:
		if m in t:
			return True
	return False

def create_test_train_lists_type(data_list):
	tellurium, selenium = [], []
	for line in data_list:
		if is_match(line,"Te"):
			tellurium.append(line)
		elif is_match(line, "Se"):
			selenium.append(line)
	return tellurium, selenium

