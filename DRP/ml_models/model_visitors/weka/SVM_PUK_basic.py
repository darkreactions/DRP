from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.AbstractModelVisitor import AbstractModelVisitor, logger
from django.core.exceptions import ImproperlyConfigured
import subprocess
import os


class SVM_PUK_basic(AbstractModelVisitor):

  maxResponseCount = 1

  def __init__(self, *args, **kwargs):
    super(SVM_PUK_basic, self).__init__(*args, **kwargs)

    self.WEKA_VERSION = "3.6" # The version of WEKA to use.
    self.PUK_OMEGA = 0.5
    self.PUK_SIGMA = 7.0

  def train(self, reactions, descriptorHeaders, filePath):
    arff_file = self._prepareArff(reactions, descriptorHeaders)

    # Currently, we support only one "response" variable.
    headers = [h for h in reactions.expandedCsvHeaders if h in descriptorHeaders]
    response_index = headers.index(list(self.statsModel.container.outcomeDescriptors)[0].csvHeader) + 1

    kernel = "\"weka.classifiers.functions.supportVector.Puk -O {} -S {}\"".format(self.PUK_OMEGA, self.PUK_SIGMA)
    command = "java weka.classifiers.functions.SMO -t {} -d {} -K {} -p 0 -c {}".format(arff_file, filePath, kernel, response_index)
    self._runWekaCommand(command)

  def predict(self, reactions, descriptorHeaders):
    arff_file = self._prepareArff(reactions, descriptorHeaders)
    model_file = self.statsModel.fileName.name

    results_file = "{}_{}.out".format(self.statsModel.pk, uuid.uuid4())
    results_path = os.path.join(settings.TMP_DIR, results_file)

    # Currently, we support only one "response" variable.
    headers = [h for h in reactions.expandedCsvHeaders if h in descriptorHeaders]
    response_index = headers.index(list(self.statsModel.container.outcomeDescriptors)[0].csvHeader) + 1

    #TODO: Validate this input.
    command = "java weka.classifiers.functions.SMO -T {} -l {} -p 0 -c {} 1> {}".format(arff_file, model_file, response_index, results_path)
    self._runWekaCommand(command)

    response = list(self.statsModel.container.outcomeDescriptors)[0]
    results = tuple((reaction, result) for reaction, result in zip(reactions, self._readWekaOutputFile(results_path)))
    return {response :results}


  def _prepareArff(self, reactions, whitelistHeaders):
    """Writes an *.arff file using the provided queryset of reactions."""
    logger.debug("Preparing ARFF file...")
    filename = "{}_{}.arff".format(self.statsModel.pk, uuid.uuid4())
    filepath = os.path.join(settings.TMP_DIR, filename)
    while os.path.isfile(filepath): #uber paranoid making sure we don't race condition
        filename = "{}_{}.arff".format(self.statsModel.pk, uuid.uuid4())
        filepath = os.path.join(settings.TMP_DIR, filename)
    with open(filepath, "w") as f:
      reactions.toArff(f, expanded=True, whitelistHeaders=whitelistHeaders)
    return filepath

  def _readWekaOutputFile(self, filename):
    """Reads a *.out file called `filename` and outputs an ordered list of the
       predicted values in that file."""
    prediction_index = 2
    with open(filename,"r") as f:
      raw_lines = f.readlines()[5:-1] # Discard the headers and ending line.
      raw_predictions = [line.split()[prediction_index] for line in raw_lines]
      predictions = [prediction.split(":")[1] for prediction in raw_predictions]
      return predictions # coerce to appropriate type later.

  def _runWekaCommand(self, command):
    """Sets the CLASSPATH necessary to use Weka, then runs a shell `command`."""
    if not settings.WEKA_PATH[self.WEKA_VERSION]:
      raise ImproperlyConfigured("'WEKA_PATH' is not set in settings.py!")

    set_path = "export CLASSPATH=$CLASSPATH:{}; ".format(settings.WEKA_PATH[self.WEKA_VERSION])
    command = set_path + command
    logger.debug("Running in Shell:\n{}".format(command))
    subprocess.check_output(command, shell=True)

