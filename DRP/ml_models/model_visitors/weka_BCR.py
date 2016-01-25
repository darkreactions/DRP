from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.AbstractModelVisitor import AbstractModelVisitor, logger
from DRP.models.descriptors import BooleanDescriptor, NumericDescriptor, CategoricalDescriptor, OrdinalDescriptor
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, OrdRxnDescriptorValue, BoolRxnDescriptorValue
from django.core.exceptions import ImproperlyConfigured
import subprocess
import os

tools=('SVM_BCR')

class SVM_BCR(AbstractModelVisitor):

  maxResponseCount = 1

  def __init__(self, *args, **kwargs):
    super(SVM_BCR, self).__init__(*args, **kwargs)

    self.WEKA_VERSION = "3.6" # The version of WEKA to use.
    self.PUK_OMEGA = 0.5
    self.PUK_SIGMA = 7.0

  def train(self, reactions, descriptorHeaders, filePath):
    arff_file = self._prepareArff(reactions, descriptorHeaders)

    # Currently, we support only one "response" variable.
    response = list(self.statsModel.container.outcomeDescriptors)[0]

    if isinstance(response, CategoricalDescriptor):
      num_classes = response.permittedValues.all().count()
      response_values = CatRxnDescriptorValue.objects().filter(reaction__in=reactions, descriptor=response)
      class_counts = [response_values.filter(value=v).count() for v in response.permittedValues().all()]
    elif isinstance(response, OrdinalDescriptor):
      num_classes = response.maximum - response.minimum + 1
      response_values = OrdRxnDescriptorValue.objects.filter(reaction__in=reactions, descriptor=response)
      class_counts = [response_values.filter(value=v).count() for v in range(response.minimum, response.maximum+1)]
    elif isinstance(response, BooleanDescriptor):
      num_classes = 2
      response_values = BoolRxnDescriptorValue.objects.filter(reaction__in=reactions, descriptor=response)
      class_counts = [response_values.filter(value=True).count(), response_values.filter(value=False).count()]
    elif isinstance(response, NumericDescriptor):
      #TODO XXX: Should this be another type of error?
      raise RuntimeError('Cannot train a classification algorithm to predict a numeric descriptor.')
    else:
      #TODO XXX: Should this be another type of error?
      raise RuntimeError('Response descriptor is not a recognized descriptor type.')
      
    # the i,j entry in cost matrix corresponds to cost of classifying an instance of class j as class i
    # since we want misclassification of an instance to be equally costly regardless of what it is classified as
    # all entries in a column not on the diagonal should be the same. The diagonal should always be 0 as correct
    # classification has no cost. So that every class is weighted the same as a whole, each class's weight is 
    # 1/(number of instances of that class). To reduce floating point arithmetic errors, we first multiply by
    # total number of data points, so each class is weighted by (total_instances)/(class_count)
    # classes for which class_count is 0 do not matter, so their weight is 0 (to avoid division by 0)
    
    total_instances = sum(class_counts)
    class_weights = [ (total_instances/float(class_count) if class_count!=0 else 0.0) for class_count in class_counts ]
    cost_matrix = [ [str(0.0) if i==j else str(class_weights[j]) for j in range(num_classes)] for i in range(num_classes) ]

    cost_matrix_string = '"[' + '; '.join([' '.join(row) for row in cost_matrix]) + ']"'
   
    # Currently, we support only one "response" variable.
    headers = [h for h in reactions.expandedCsvHeaders if h in descriptorHeaders]
    response_index = headers.index(list(self.statsModel.container.outcomeDescriptors)[0].csvHeader) + 1

    kernel = '"weka.classifiers.functions.supportVector.Puk -O {} -S {}"'.format(self.PUK_OMEGA, self.PUK_SIGMA)
    command = "java weka.classifiers.meta.CostSensitiveClassifier -cost-matrix {} -W weka.classifiers.functions.SMO -t {} -d {} -p 0 -c {} -- -K {}".format(cost_matrix_string, arff_file, filePath, response_index, kernel)

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

