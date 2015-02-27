from django.db import models
from DRP.settings import MODEL_DIR

class ModelStats(models.Model):
  class Meta:
    app_label = "DRP"

  # Model Statistics
  confusion_table = models.TextField(default="{}")
  train_confusion_table = models.TextField(default="{}")
  headers = models.TextField(default="[]")
  correct_vals = models.CharField("Correct Values", max_length=100,
                                  default="[\"3\",\"4\"]")

  # Model Descriptors
  title = models.CharField("Title", max_length=100, default="untitled")
  description = models.TextField(default="")
  tags = models.TextField(default="")
  iterations = models.IntegerField("iterations", default=1)

  # Model Status and Location
  filename = models.CharField("Filename", max_length=128,
                                          default=MODEL_DIR+"untitled.model")
  active = models.BooleanField("Active", default=True)
  start_time = models.DateTimeField(blank=True, null=True)
  end_time = models.DateTimeField(blank=True, null=True)
  usable = models.BooleanField("Usable", default=True)

  # Available model types.
  library = models.CharField("Library", max_length=128, default="weka")
  tool = models.CharField("Tool", max_length=128, default="svc")
  response = models.CharField("Response", max_length=128, default="outcome")



  def construct(self, title, data, description="", tags="",
                                     library="weka", tool="svc",
                                     response="outcome",
                                     filename="", force=False,
                                     usable=True, active=False,
                                     preprocessor=None, postprocessor=None,
                                     splitter=None,
                                     test_set=None,
                                     clean_tmp_files=True,
                                     debug=False):
    """
    Using the headers and rows specified in `data`, train and test a model.

    Description of the options are below:

    force -- Forces the overwrite of any models using the same filepath.
    usable -- Whether this model can be used for predictions.
    active -- Whether this model is viewable in the dashboard.
    debug -- Enables/disables helpful stdout messages. Does not affect errors.
    clean_tmp_files -- Enable/disable the auto-deletion of temporary files
                       that are generated (Eg: *.arff, *.out).

    preprocessor -- A function that will be passed the uncleaned data
                    *before* splitting.
    postprocessor -- A function that will be passed the data returned from
                     the `splitter` function and the headers. Should return
                     a tuple with (data, headers).
    splitter -- A function that may be passed the kwarg "headers" and the
                preprocessed data. Should return a dictionary with the split
                data (eg: {"all":data_1, "train":data_2, "test":data_3}).
    """

    from DRP.fileFunctions import file_exists
    import datetime, random

    # Use a custom splitter-function if specified.
    if not splitter:
      from DRP.model_building.splitters import default_splitter as splitter

    # Make sure we can actually write the model file.
    if not filename:
      filename = "".join(filter(str.isalnum, title))
    filename += "__{}__{}".format(library, tool)
    self.filename = filename

    # Don't overwrite models unless "force=True" is specified.
    if file_exists(self.get_path()) and not force:
      message = "Model '{}' already exists: use 'force=True'".format(self.get_path())
      raise Exception(message)

    if file_exists(self.get_path()) and force:
      if debug:
        print "Forcing file overwrite: {}".format(self.get_path())

      # Disable any models that point to this file since they can no longer be used.
      others = get_models_from_filename(self.filename)
      others.update(usable = False)


    self.start_time = datetime.datetime.now()

    # Save the description and fields of this model.
    self.title = title
    self.description = description
    self.tags = tags
    self.usable = usable
    self.active = active
    self.library = library
    self.response = response
    self.tool = tool
    self.filename = filename

    # Don't let the original data be modified in any way.
    data = data[:]

    if debug:
      print "Starting generation of '{}' using '{}' on {} entries".format(self.tool, self.library, len(data))

    # Pre-process the data if it is not already processed.
    if preprocessor and not test_set:
      if debug: print "Pre-processing... ({})".format(preprocessor.__name__)
      data = preprocessor(data)

    headers = data.pop(0)
    random.shuffle(data) #Randomize the data before splitting.

    if not test_set:
      # Split the data.
      if debug: print "Splitting data... ({})".format(splitter.__name__)
      split_data = splitter(data, headers=headers)
      if debug:
        splits = {key:len(val) for key, val in split_data.items()}
        print "Splits: {}".format(splits)
    else:
      if debug: print "Using predefined test data..."
      split_data = {"all":data+test_set, "test":test_set, "train":data}


    # If specified, post-process the data after splitting.
    if postprocessor:
      if debug: print "Post-processing... ({})".format(postprocessor.__name__)
      split_data, headers = postprocessor(split_data, headers)

    # Save the headers now that the data has successfully been split.
    if self.response not in headers:
      raise Exception("Response '{}' not found in headers!".format(self.response))

    self.set_headers(headers)

    if debug:
      print "Using {} headers...".format(len(headers))

    # Get the temporary value-sets of each field.
    self.val_map = self._get_val_map(split_data["all"])

    # Train/fit the model
    if debug: print "Training model..."
    self._train_model(split_data["train"], debug=debug)


    # Set the confusion table for a given data-set.
    if debug: print "Testing model..."
    self._test_model(split_data["test"], debug=debug, table="test")

    self.end_time = datetime.datetime.now()
    self.save()

    if clean_tmp_files:
      self._delete_tmp_files()

    if debug:
      print "Complete! Took {} seconds.".format(self._construction_time())

    return self


  def _construction_time(self):
    return self.end_time-self.start_time


  def _get_val_map(self, data):
    fields = self.get_headers()
    val_dict = {field:{row[i] for row in data} for i, field in enumerate(fields)}
    return val_dict

  def get_path(self):
    from DRP.settings import MODEL_DIR
    from DRP.fileFunctions import createDirIfNecessary

    createDirIfNecessary(MODEL_DIR)

    return MODEL_DIR + self.filename

  def get_headers(self):
    import json
    return json.loads(self.headers)



  def _strip_response(self, data):
    """
    Separates the predictors in a data matrix from the response.
    """

    if len(data)==0:
      raise Exception("`data` cannot be empty!")

    headers = self.get_headers()
    resp = headers.index(self.response)

    preds = [[elem for i, elem in enumerate(row) if i!=resp] for row in data]
    resps = [row[resp] for row in data]

    return preds, resps

  def _read_weka_output(self, filepath):
    def _is_number(elem):
      try:
        float(elem)
        return True
      except:
        return False

    predicted_index = 2
    with open(filepath) as f:
      # Get the `raw_predictions` from the file, but ignore the headings/footings.
      content = f.readlines()
      raw_predictions = content[5:-1]

      predictions = [entry.split()[predicted_index] for entry in raw_predictions]

      # WEKA's predictions are in the format: "class_num:value".
      # Note that class_num is arbitrary, thus we only need "value".
      predictions = [prediction.split(":")[1] for prediction in predictions]

      return predictions


  def predict(self, predictors, debug=False, table="test"):
    from DRP.fileFunctions import get_django_path, file_exists
    import subprocess

    if not self.usable:
      raise Exception("ModelStats is set to 'unusable' and cannot use `predict`.")

    if not file_exists(self.get_path()):
      raise Exception("Filepath to model is invalid: '{}'".format(self.get_path()))


    if self.library == "weka":
      arff_path = self._make_arff(table, predictors, debug=debug)

      # Results path will be named *.out instead of *.arff.
      results_path = "out".join(arff_path.rsplit("arff", 1))

      move = "cd {};".format(get_django_path())
      comm = "bash DRP/model_building/make_predictions.sh"
      args = " {} {} {}".format(arff_path, self.get_path(), results_path)

      subprocess.check_output(move+comm+args, shell=True)
      predictions = self._read_weka_output(results_path)

    elif self.library == "sklearn":
      model = self._load_model()
      predictions = model.predict(predictors)

    else:
      raise Exception("Unknown library specified in predict!")


    if self.library in {"sklearn", "weka"}:
      predictions = map(int, map(float, predictions))

    return predictions


  def _load_model(self):
    from sklearn.externals import joblib

    if self.library == "sklearn":
      return joblib.load(self.get_path())

    else:
      raise Exception("Illegal model library specified! Aborting file-load!")


  def _test_model(self, data, debug=False, table="test"):
    from DRP.model_building.confusion_table import make_confusion_dict

    if self.library=="weka":
      predictors = data
      index = self.get_headers().index(self.response)
      responses = [row[index] for row in data]

    else:
      predictors, responses = self._strip_response(data)

    guesses = self.predict(predictors, debug=debug, table=table)

    if self.library in {"sklearn", "weka"}:
      responses = map(int, map(float, responses))

    cm = make_confusion_dict(guesses, responses)
    self.set_confusion_table(cm, table=table)



  def _train_model(self, data, debug=False):
    from DRP.model_building.sklearn_methods import get_model
    from DRP.fileFunctions import get_django_path
    from sklearn.externals import joblib
    import subprocess

    if self.library=="weka":
      if self.tool!="svc":
        raise Exception("Only svc supported in WEKA currently.")

      arff_name = self._make_arff("train", data, debug=debug)

      move = "cd {}; ".format(get_django_path())
      command = "bash DRP/model_building/make_model.sh "
      args = "{} {}".format(self.get_path(), arff_name)

      if debug:
        print "\tShell: {} {} {}".format(move, command, args)

      output = subprocess.check_output(move + command + args, shell=True,
                                       stderr=subprocess.STDOUT)

      # If WEKA raised an error, throw that error.
      if "Exception" in output:
        raise Exception(output)

      if debug:
        print "Testing against training data..."
      self._test_model(data, debug=debug, table="train")

    elif self.library=="sklearn":
      model, description = get_model(self.tool)
      predictors, responses = self._strip_response(data)
      model.fit(predictors, responses)

      # Save the sklearn model to a file (with maximum compression).
      joblib.dump(model, self.get_path(), compress=9)

    else:
      raise Exception("_train_model called with illegal library!")

  def _get_weka_value_map(self):
    def _is_string(elem):
      return type(elem)==str or type(elem)==unicode

    val_dict = {}

    for field, val_set in self.val_map.items():
      if any(map(_is_string, val_set)):
        if field in self.val_map:
          val_set = self.val_map[field]

        # EG: convert {"hello", 1, u"world"} to {hello,1,world}.
        innards = ",".join(map(str,val_set)).replace("\"","")
        val_dict[field] = "{"+innards+"}"

      else:
        val_dict[field] = "NUMERIC"

    return val_dict

  def _weka_header(self, data):
    headers = self.get_headers()
    value_map = self._get_weka_value_map()

    res = "%  COMMENT \n%  NAME, DATE\n@relation rec_system"

    for header in headers:
      res += "\n@ATTRIBUTE " + header + " " + value_map[header]

    res += "\n\n@DATA\n"
    return res

  def _delete_tmp_files(self):
    from DRP.settings import TMP_DIR
    from DRP.fileFunctions import delete_tmp_file

    # Variable Setup
    auto_deleted_keys = ["test", "train"]
    auto_deleted_exts = ["arff", "out"]

    for ext in auto_deleted_exts:
      for key in auto_deleted_keys:
        filepath = "{}{}_{}.{}".format(TMP_DIR, self.filename, key, ext)
        delete_tmp_file(filepath)


  def _make_arff(self, key, data, debug=False):
    """
    Saves the data to an ARFF file.

    Numeric values are given NUMERIC headings and non-numeric values
    are treated as classes. Passing a value as a `str` will cause it to
    be treated as a class (such as for "outcome").
    """

    from DRP.fileFunctions import createDirIfNecessary
    from DRP.settings import TMP_DIR

    createDirIfNecessary(TMP_DIR)
    full_name = "{}{}_{}.arff".format(TMP_DIR, self.filename, key)

    with open(full_name, "w") as f:
      header_content = self._weka_header(data)
      f.write(header_content + "\n")
      for row in data:
        f.write(",".join([str(entry) for entry in row]) + "\n")

    if debug: print "\tARFF: {}".format(full_name)
    return full_name


  def set_correct_vals(self, correct_list):
    import json
    if not correct_list:
      correct_list = ["4","3"]
    self.correct_vals = json.dumps(correct_list)
    self.save()


  def load_all_vals(self):
    return sorted(self.load_confusion_dict().keys())

  def load_correct_vals(self):
    import json
    return sorted(json.loads(self.correct_vals))

  def load_incorrect_vals(self):
    correct = set(self.load_correct_vals())
    all_vals = self.load_all_vals()
    return [val for val in all_vals if val not in correct]


  def check_usability(self):
    import os

    try:
      assert( len(self.load_all_vals()) > 0 )
      assert( len(self.load_correct_vals()) > 0 )
      assert( os.path.isfile(self.get_path()) )
      usable = True
    except:
      usable = False

    if usable!= self.usable:
      self.usable = usable
      self.save()

    return usable


  def set_headers(self, field_list):
    import json
    self.headers = json.dumps(field_list)
    self.save()


  def set_confusion_table(self, conf_json, table="test"):
    import json

    cm = json.dumps(conf_json)

    if table=="test":
      self.confusion_table = cm
    elif table=="train":
      self.train_confusion_table = cm
    else:
      raise Exception("Illegal confusion table specified!")

    self.save()

  def graph_confusion_table(self, table="test"):
    import matplotlib.pyplot as plt
    ticks = self.load_all_vals()

    raw_cm = self.load_confusion_table(normalize=False, headers=False, table=table)
    cm = [map(float,row) for row in raw_cm]

    plt.matshow(cm, cmap=plt.cm.OrRd)
    plt.colorbar()
    plt.ylabel('True')
    plt.xlabel('Predicted')
    plt.xticks(range(len(ticks)), ticks)
    plt.yticks(range(len(ticks)), ticks)

    if self.title:
      plt.title(self.title)

    plt.show()


  def load_confusion_dict(self, table="test"):
    import json

    if table=="test":
      cm = self.confusion_table
    elif table=="train":
      cm = self.train_confusion_table
    else:
      raise Exception("Illegal confusion table specified!")

    return json.loads(cm)

  def load_confusion_table(self, normalize=True, headers=True, table="test"):
    """
    Confusion Dict:
    Abstract Format:
    {
      "Actual Value": {"Predicted Val":amount, "Predicted Val2": amount2, ... }
      ...
    }

    Actual Format:
    {
      "1":{"1":#, "2":#, "3":#, "4":#},
      "2":{"1":#, "2":#, "3":#, "4":#},
      "3":{"1":#, "2":#, "3":#, "4":#},
      "4":{"1":#, "2":#, "3":#, "4":#},
    }

    Confusion Matrix:
            __Predicted___
      __    V1   V2   V3 ...
      A  V1 #    #    #
      c
      t  V2 #    #    #
      u
      a  V3 #    #    #
      l
      __ ...
    """

    try:
      confusion_dict = self.load_confusion_dict(table=table)
    except:
      return []

    values = sorted(confusion_dict.keys())

    matrix = [[""]+values] if headers else []

    denom = self.total() if normalize else 1

    for value in values:
      guess_dict = confusion_dict[value]
      row = [value] if headers else []
      row += [guess_dict[v]/denom if v in guess_dict else 0 for v in values]
      matrix.append(row)

    return matrix


  def count(self, normalize=False, guesses=None, actuals=None, ranges=True, false_guess=False, table="test"):
    # Variable Setup
    conf_table = self.load_confusion_table(normalize=normalize, table=table)

    guess_headers = conf_table.pop(0)[1:] # Remove the empty cell in [0,0].
    actual_headers = [row.pop(0) for row in conf_table]

    if not actuals: actuals = actual_headers
    if not guesses: guesses = guess_headers

    c = 0

    for value in (actuals+guesses):
      if value not in guess_headers or value not in actual_headers:
        raise Exception("Value '{}' not found in confusion table!".format(value))

    for i, guess in enumerate(guess_headers):
      for j, actual in enumerate(actual_headers):
        if ranges:
          if false_guess:
            if actual in actuals and guess in guesses and actual not in guesses:
              c += conf_table[j][i]

          else:
            if actual in actuals and guess in guesses:
              c += conf_table[j][i]

        else:
          if false_guess:
            if actual!=guess and actual in actuals and guess in guesses:
              c += conf_table[j][i]

          else:
            if actual==guess and actual in actuals and guess in guesses:
              c += conf_table[j][i]

    return c


  def total(self, table="test"):
    conf_table = self.load_confusion_dict(table=table)
    int_total = sum([int(val) for correct,guesses in conf_table.items()
                              for key,val in guesses.items()])
    return float(int_total)

  # Convenience Wrappers
  def true_positives(self, ranges=True, normalize=False, table="test"):
    corrects = self.load_correct_vals()
    return self.count(guesses=corrects, actuals=corrects,
                      normalize=normalize, ranges=ranges, table=table)

  def true_negatives(self, ranges=True, normalize=False, table="test"):
    incorrects = self.load_incorrect_vals()
    return self.count(guesses=incorrects, actuals=incorrects,
                      normalize=normalize, ranges=ranges, table=table)

  def false_positives(self, ranges=True, normalize=False, table="test"):
    corrects = self.load_correct_vals()
    return self.count(guesses=corrects, false_guess=True,
                      normalize=normalize, ranges=ranges, table=table)

  def false_negatives(self, ranges=True, normalize=False, table="test"):
    incorrects = self.load_incorrect_vals()
    return self.count(guesses=incorrects, false_guess=True,
                      normalize=normalize, ranges=ranges, table=table)


  def accuracy(self, ranges=True, table="test"):
    denom = float(self.total(table=table))
    if denom:
      tp = self.true_positives(ranges=ranges,table=table)
      tn = self.true_negatives(ranges=ranges,table=table)
      return (tp + tn)/denom
    else:
      return 0

  def precision(self, ranges=True, table="test"):
    tp = self.true_positives(ranges=ranges,table=table)
    fp = self.false_positives(ranges=ranges,table=table)
    denom = float(tp + fp)
    if denom:
      return tp/denom
    else:
      return 0

  def recall(self, ranges=True, table="test"):
    tp = self.true_positives(ranges=ranges,table=table)
    fn = self.false_negatives(ranges=ranges,table=table)
    denom = float(tp + fn)
    if denom:
      return tp/denom
    else:
      return 0

  def user_satisfaction(self):
    from DRP.models import Recommendation
    recs = Recommendation.objects.filter(model_version=self)
    if recs.exists():
      return recs.filter(nonsense=False).count()/float(recs.count())
    else:
      return 0

  def pvalue(self):
    #TODO: return the p-value for this model.
    return float("inf") #TODO: Not this.

  def pretty_stats(self, classes=None):
    print "Model: {}".format(self.title)
    print "Description: {}".format(self.description)
    raw = self.stats()
    for categorization, stat_dict in raw.items():
      print "{}-categorization:".format(categorization)
      for stat, val in stat_dict.items():
        print "\t{}: {}".format(stat, val)


  def stats(self, category=None):

    tests = {
      "2-test": {
        "Test Size":self.total(table="test"),
        "Train Size":self.total(table="train"),
        "Accuracy":self.accuracy(ranges=True),
        "Precision":self.precision(ranges=True),
        "Recall":self.recall(ranges=True),
        "% TP":self.true_positives(normalize=True, ranges=True),
        "% FP":self.false_positives(normalize=True, ranges=True),
        "% TN":self.true_negatives(normalize=True, ranges=True),
        "% FN":self.false_negatives(normalize=True, ranges=True),
      },
      "4-test": {
        "Test Size":self.total(table="test"),
        "Train Size":self.total(table="train"),
        "Accuracy":self.accuracy(ranges=False),
        "Precision":self.precision(ranges=False),
        "Recall":self.recall(ranges=False),
        "% TP":self.true_positives(normalize=True, ranges=False),
        "% FP":self.false_positives(normalize=True, ranges=False),
        "% TN":self.true_negatives(normalize=True, ranges=False),
        "% FN":self.false_negatives(normalize=True, ranges=False),
      },
      "2-train": {
        "Test Size":self.total(table="test"),
        "Train Size":self.total(table="train"),
        "Accuracy":self.accuracy(ranges=True, table="train"),
        "Precision":self.precision(ranges =True, table="train"),
        "Recall":self.recall(ranges=True, table="train"),
        "% TP":self.true_positives(normalize=True, ranges=True, table="train"),
        "% FP":self.false_positives(normalize=True, ranges=True, table="train"),
        "% TN":self.true_negatives(normalize=True, ranges=True, table="train"),
        "% FN":self.false_negatives(normalize=True, ranges=True, table="train"),
      },
      "4-train": {
        "Test Size":self.total(table="test"),
        "Train Size":self.total(table="train"),
        "Accuracy":self.accuracy(ranges=False, table="train"),
        "Precision":self.precision(ranges=False, table="train"),
        "Recall":self.recall(ranges=False, table="train"),
        "% TP":self.true_positives(normalize=True, ranges=False, table="train"),
        "% FP":self.false_positives(normalize=True, ranges=False, table="train"),
        "% TN":self.true_negatives(normalize=True, ranges=False, table="train"),
        "% FN":self.false_negatives(normalize=True, ranges=False, table="train"),
      },
    }

    if category:
      if type(category) in {set, list}:
        tests = {key:val for key,val in tests if key in category}

      else:
        tests = tests[category]

    return tests


  def print_confusion_table(self, normalize=True, table="test"):
    def truncate_floats(row):
      cleaned = []
      for elem in row:
        try:
          cleaned.append("{:.3f}".format(elem))
        except:
          cleaned.append(elem)
      return cleaned

    print "{} Confusion Table:".format(table.capitalize())
    conf_table = self.load_confusion_table(normalize=normalize, table=table)
    if conf_table:
      heading = "(%)" if normalize else ""
      print "\t\t\tPredicted {}".format(heading)
      for row in conf_table:
        cleaned_row = map(str, truncate_floats(row) )
        print "\t"+"\t".join( cleaned_row )
    else:
      print "\t[ Confusion Matrix Unavailable ]"


  def print_model_info(self, prefix="\t"):
    print prefix+"Name: '{}'".format(self.title)
    print prefix+"Description: '{}'".format(self.description)
    print prefix+"ID: {}".format(self.id)
    print prefix+"Finished: {}".format(self.end_time)
    print prefix+"Construction Time: {}".format(self._construction_time())
    print
    print prefix+"Filename: '{}'".format(self.filename)
    print prefix+"Library: '{}'".format(self.library)
    print prefix+"Tool: '{}'".format(self.tool)
    print prefix+"# Headers: '{}'".format(len(self.get_headers()))
    print prefix+"# Iterations: '{}'".format(self.iterations)
    print prefix+"Correct Values: {}".format(self.load_correct_vals())

  def summary(self, pre="\t"):
    self.print_model_info()
    print ""
    self.print_confusion_table(table="train")
    print ""
    self.print_confusion_table(table="test")
    print ""

    stats = self.stats().items()
    stats.sort(key=lambda tup: tup[0])
    for stat, key in stats:
      print "{}{}: {}".format(pre, stat, key)


  def __unicode__(self):
    return "Model {} (Active:{}; Usable:{})".format(self.title,
                                                    self.active,
                                                    self.usable)

  def __repr__(self):
    return self.__unicode__()


def get_models_from_filename(filename):
  return ModelStats.objects.filter(filename=filename)


