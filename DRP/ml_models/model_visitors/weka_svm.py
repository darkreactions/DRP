from django.conf import settings
from AbstractModelVisitor import AbstractModelVisitor
from django.core.exceptions import ImproperlyConfigured
import subprocess
import time
import os
import datetime

class ModelVisitor(AbstractModelVisitor):

  def __init__(self, modelContainer, **kwargs):
    super(ModelVisitor, self).__init__("weka", "svm", 1, modelContainer, **kwargs)

    self.model_filename = "{}.model".format(self.getModelTag())
    self.WEKA_VERSION = "3.6" # The version of WEKA to use.
    self.debug = False # Set to "True" to enable printing of debug-messages.

  def _train(self):
    reactions = self.getTrainingData()
    model_filepath = os.path.join(settings.MODEL_DIR, self.model_filename)
    arff_file = self._prepareArff(reactions, suffix="train")

    self.stats_model.start_time = datetime.datetime.now()

    kernel = "\"weka.classifiers.functions.supportVector.Puk -O 0.5 -S 7\""
    command = "java weka.classifiers.functions.SMO -t {} -d {} -K {} -p 0".format(arff_file, model_filepath, kernel)
    self._runWekaCommand(command)

    self.setModelFile(model_filepath)

    self.stats_model.end_time = datetime.datetime.now()
    self.stats_model.save()

  def predict(self, reactions, suffix="predict"):
    arff_file = self._prepareArff(reactions, suffix=suffix)
    model_file = self.getModelFilename()

    results_file =  "{}_{}.out".format(self.getModelTag(), suffix)
    results_path = os.path.join(settings.TMP_DIR, results_file)

    # Currently, we support only one "response" variable.
    response_index = self.getPredictors().count()+1

    #TODO: Validate this input.
    command = "java weka.classifiers.functions.SMO -T {} -l {} -p 0 -c {} 1> {}".format(arff_file, model_file, response_index, results_path)
    self._runWekaCommand(command)

    response = self.getResponses()[0]
    return { response : self._readWekaOutputFile(results_path) }


  def _prepareArff(self, reactions, suffix=""):
    """Writes an *.arff file using the provided queryset of reactions."""
    if self.debug: print "Preparing ARFF file..."

    filename = "{}_{}_{}.arff".format(self.getModelTag(), suffix, time.time())
    filepath = os.path.join(settings.TMP_DIR, filename)
    with open(filepath, "w") as f:
      whitelist = list(self.getPredictors()) + list(self.getResponses())
      reactions.toArff(f, expanded=True, whitelistDescriptors=whitelist)

    return filepath

  def _readWekaOutputFile(self, filename):
    """Reads a *.out file called `filename` and outputs an ordered list of the
       predicted values in that file."""
    prediction_index = 2
    with open(filename,"r") as f:
      raw_lines = f.readlines()[5:-1] # Discard the headers and ending line.
      raw_predictions = [line.split()[prediction_index] for line in raw_lines]
      predictions = [prediction.split(":")[1] for prediction in raw_predictions]
      return map(float, predictions)

  def _runWekaCommand(self, command):
    """Sets the CLASSPATH necessary to use Weka, then runs a shell `command`."""
    if not settings.WEKA_PATH[self.WEKA_VERSION]:
      raise ImproperlyConfigured("'WEKA_PATH' is not set in settings.py!")

    set_path = "export CLASSPATH=$CLASSPATH:{}; ".format(settings.WEKA_PATH[self.WEKA_VERSION])
    command = set_path + command

    if self.debug: print "Running in Shell:\n{}".format(command)

    subprocess.check_output(command, shell=True)
