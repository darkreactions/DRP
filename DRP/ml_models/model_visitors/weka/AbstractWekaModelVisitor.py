from django.conf import settings
import uuid
from DRP.models import rxnDescriptors
from DRP.ml_models.model_visitors.AbstractModelVisitor import AbstractModelVisitor, logger
from DRP.models.descriptors import BooleanDescriptor, NumericDescriptor, CategoricalDescriptor, OrdinalDescriptor
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, OrdRxnDescriptorValue, BoolRxnDescriptorValue
from django.core.exceptions import ImproperlyConfigured
import subprocess
import os
from abc import abstractmethod, abstractproperty
import warnings
from itertools import chain

class AbstractWekaModelVisitor(AbstractModelVisitor):

    maxResponseCount = 1
    WEKA_VERSION = "3.6" # The version of WEKA to use.

    def __init__(self, BCR=False, *args, **kwargs):
        self.BCR = BCR
        
        super(AbstractWekaModelVisitor, self).__init__(*args, **kwargs)

        # This is a bit hackier, but I don't think anything like abstractattribute is implemented in abc
        try:
            self.wekaCommand
        except AttributeError:
            raise NotImplementedError('Subclasses of AbstractWekaModelVisitor must define wekaCommand')

    def _prepareArff(self, reactions, whitelistHeaders, verbose=False):
        """Writes an *.arff file using the provided queryset of reactions."""
        logger.debug("Preparing ARFF file...")
        filename = "{}_{}.arff".format(self.statsModel.pk, uuid.uuid4())
        filepath = os.path.join(settings.TMP_DIR, filename)
        while os.path.isfile(filepath):  # uber paranoid making sure we don't race condition
            filename = "{}_{}.arff".format(self.statsModel.pk, uuid.uuid4())
            filepath = os.path.join(settings.TMP_DIR, filename)
        if verbose:
            print "Writing arff to {}".format(filepath)
        with open(filepath, "w") as f:
            reactions.toArff(f, expanded=True, whitelistHeaders=whitelistHeaders)
        return filepath

    def _readWekaOutputFile(self, filename, typeConversionFunction):
        """
        Reads a *.out file called `filename` and outputs an ordered list of the
        predicted values in that file.
        """
        prediction_index = 2
        with open(filename, "r") as f:
            raw_lines = f.readlines()[5:-1]  # Discard the headers and ending line.
            raw_predictions = [line.split()[prediction_index] for line in raw_lines]
            predictions = [typeConversionFunction(prediction) for prediction in raw_predictions]
        return predictions

    def _runWekaCommand(self, command, verbose=False):
        """Sets the CLASSPATH necessary to use Weka, then runs a shell `command`."""
        if not settings.WEKA_PATH[self.WEKA_VERSION]:
            raise ImproperlyConfigured("'WEKA_PATH' is not set in settings.py!")
        set_path = "export CLASSPATH=$CLASSPATH:{}; ".format(settings.WEKA_PATH[self.WEKA_VERSION])
        command = set_path + command
        logger.debug("Running in Shell:\n{}".format(command))
        if verbose:
            print "Running in Shell:\n{}".format(command)
        subprocess.check_output(command, shell=True)
        # TODO XXX. figure out some way to throw an error if weka errors (sends stuff to stderr. weka does not generate proper return codes)
        # This broke things in some cases for reasons not clear to me
        # wekaProc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True, shell=True) 
        # wekaProc.wait()
        # if wekaProc.returncode == 0:
        #     res, resErr = wekaProc.communicate()
        #     if resErr:
        #         raise RuntimeError("Weka returned an error: {}".format(resErr))
        # else:
        #     raise RuntimeError("Weka returned nonzero exit code")


    def BCR_cost_matrix(self, reactions, response):
        if isinstance(response, CategoricalDescriptor):
            num_classes = response.permittedValues.all().count()
            response_values = CatRxnDescriptorValue.objects().filter(reaction__in=reactions, descriptor=response)
            # TODO XXX Make this happen in the right order
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
            
            
        ## the i,j entry in cost matrix corresponds to cost of classifying an instance of class i as class j
        ## since we want misclassification of an instance to be equally costly regardless of what it is classified as
        ## all entries in a row not on the diagonal should be the same. The diagonal should always be 0 as correct
        ## classification has no cost. So that every class is weighted the same as a whole, each class's weight is 
        ## 1/(number of instances of that class). To reduce floating point arithmetic errors, we first multiply by
        ## total number of data points, so each class is weighted by (total_instances)/(class_count)
        ## classes for which class_count is 0 do not matter, so their weight is 0 (to avoid division by 0)
        ## For boolean classification, True is class 0 and False is class 1 (Because that's how it's set up in the toArff function)
        ## TODO XXX Actually check the order of classes from the arff file?
        
        total_instances = sum(class_counts)
        class_weights = [ (total_instances/float(class_count) if class_count!=0 else 0.0) for class_count in class_counts ]
        cost_matrix = [ [str(0.0) if i==j else str(class_weights[i]) for j in range(num_classes)] for i in range(num_classes) ]

        cost_matrix_string = '"[' + '; '.join([' '.join(row) for row in cost_matrix]) + ']"'

        return cost_matrix_string

    def wekaTrainOptions(self):
        """
        Returns any additional commands specific to the classifier
        """
        return ""

    def train(self, verbose=False):
        reactions = self.statsModel.trainingSet.reactions.all()
        # TODO XXX update this to make use of csvHeader as a query
        descriptorHeaders = [d.csvHeader for d in chain(self.statsModel.container.descriptors, self.statsModel.container.outcomeDescriptors)]
        filePath = self.statsModel.outputFile.name
        if not self.statsModel.inputFile.name:
            self.statsModel.inputFile = self._prepareArff(reactions, descriptorHeaders, verbose)
            self.statsModel.save(update_fields=['inputFile'])
        elif not os.path.isfile(self.statsModel.inputFile.name):
            if self.invalid:
                raise RuntimeError('Could not find statsModel arff file and model is invalid')
            else:
                raise warning.warn('Could not find statsModel arff file, but model is valid, so recreating')
                self.statsModel.inputFile.name = self._prepareArff(reactions, descriptorHeaders, verbose)
                self.statsModel.save(update_fields=['inputFile'])

        arff_file = self.statsModel.inputFile.name

        # Currently, we support only one "response" variable.
        response = list(self.statsModel.container.outcomeDescriptors)[0]
        headers = [h for h in reactions.expandedCsvHeaders() if h in descriptorHeaders]
        response_index = headers.index(response.csvHeader) + 1

        if self.BCR:
            cost_matrix_string = self.BCR_cost_matrix(reactions, response)
            command = "java weka.classifiers.meta.CostSensitiveClassifier -cost-matrix {} -W {} -t {} -d {} -p 0 -c {} -- {}".format(cost_matrix_string, self.wekaCommand, arff_file, filePath, response_index, self.wekaTrainOptions())
        else:
            command = "java {} -t {} -d {} -p 0 -c {} {}".format(self.wekaCommand, arff_file, filePath, response_index, self.wekaTrainOptions())
        self._runWekaCommand(command, verbose=verbose)

    def predict(self, reactions, verbose=False):
        descriptorHeaders = [d.csvHeader for d in chain(self.statsModel.container.descriptors, self.statsModel.container.outcomeDescriptors)]
        
        arff_file = self._prepareArff(reactions, descriptorHeaders, verbose=verbose)
        model_file = self.statsModel.outputFile.name

        results_file = "{}_{}.out".format(self.statsModel.pk, uuid.uuid4())
        results_path = os.path.join(settings.TMP_DIR, results_file)

        # Currently, we support only one "response" variable.
        headers = [h for h in reactions.expandedCsvHeaders() if h in descriptorHeaders]
        response = list(self.statsModel.container.outcomeDescriptors)[0]
        response_index = headers.index(response.csvHeader) + 1

        # TODO: Validate this input.
        command = "java {} -T {} -l {} -p 0 -c {} 1> {}".format(self.wekaCommand, arff_file, model_file, response_index, results_path)
        if verbose:
            print "Writing results to {}".format(results_path)
        self._runWekaCommand(command, verbose=verbose)

        if isinstance(response, rxnDescriptors.BoolRxnDescriptor):
            typeConversionFunction = booleanConversion
        elif isinstance(response, rxnDescriptors.OrdRxnDescriptor):
            typeConversionFunction = ordConversion
        elif isinstance(response, rxnDescriptors.NumRxnDescriptor):
            typeConversionFunction = numConversion
        elif isinstance(response, rxnDescriptors.CatRxnDescriptor):
            typeConversionFunction = str
        else:
            raise TypeError("Response descriptor is of invalid type {}".format(type(response)))
        results = tuple((reaction, result) for reaction, result in zip(reactions, self._readWekaOutputFile(results_path, typeConversionFunction)))
        return {response: results}


def numConversion(s):
    return float(s)

def ordConversion(s):
    s = s.split(':')[1]
    return int(s)

def booleanConversion(s):
    s = s.split(':')[1]

    if s.lower() == 'true':
        return True
    elif s.lower() == 'false':
        return False
    else:
        raise ValueError("Tried to convert string to boolean when string was neither 'True' nor 'False' but {}".format(s))
