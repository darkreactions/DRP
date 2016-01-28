from django.conf import settings
import uuid
from DRP.ml_models.model_visitors.weka.AbstractWekaModelVisitor import AbstractWekaModelVisitor
from DRP.models.descriptors import BooleanDescriptor, NumericDescriptor, CategoricalDescriptor, OrdinalDescriptor
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, OrdRxnDescriptorValue, BoolRxnDescriptorValue
import os


class SVM_PUK_BCR(AbstractWekaModelVisitor):

  maxResponseCount = 1

  def __init__(self, *args, **kwargs):
    super(SVM_PUK_BCR, self).__init__(*args, **kwargs)

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
      raise TypeError('Cannot train a classification algorithm to predict a numeric descriptor.')
    else:
      raise TypeError('Response descriptor is not a recognized descriptor type.')
      
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
