from django.conf import settings
import uuid
from DRP.models import rxnDescriptors
from DRP.ml_models.model_visitors.AbstractModelVisitor import AbstractModelVisitor, logger
from django.core.exceptions import ImproperlyConfigured
import subprocess
import os
from abc import abstractmethod


class AbstractWekaModelVisitor(AbstractModelVisitor):

    maxResponseCount = 1

    def __init__(self, *args, **kwargs):
        super(AbstractWekaModelVisitor, self).__init__(*args, **kwargs)

        self.WEKA_VERSION = "3.6"  # The version of WEKA to use.

    def _prepareArff(self, reactions, whitelistHeaders):
        """Writes an *.arff file using the provided queryset of reactions."""
        logger.debug("Preparing ARFF file...")
        filename = "{}_{}.arff".format(self.statsModel.pk, uuid.uuid4())
        filepath = os.path.join(settings.TMP_DIR, filename)
        while os.path.isfile(filepath):  # uber paranoid making sure we don't race condition
            filename = "{}_{}.arff".format(self.statsModel.pk, uuid.uuid4())
            filepath = os.path.join(settings.TMP_DIR, filename)
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
            predictions = [typeConversionFunction(prediction.split(":")[1]) for prediction in raw_predictions]
        return predictions

    def _runWekaCommand(self, command):
        """Sets the CLASSPATH necessary to use Weka, then runs a shell `command`."""
        if not settings.WEKA_PATH[self.WEKA_VERSION]:
            raise ImproperlyConfigured("'WEKA_PATH' is not set in settings.py!")
        set_path = "export CLASSPATH=$CLASSPATH:{}; ".format(settings.WEKA_PATH[self.WEKA_VERSION])
        command = set_path + command
        logger.debug("Running in Shell:\n{}".format(command))
        subprocess.check_output(command, shell=True)

    @abstractmethod
    def wekaPredict(self, arff_file, model_file, response_index, results_path):
        """
        A function meant to be overridden by actual ModelVisitor classes.
        Run the appropriate weka prediction command.
        """

    def predict(self, reactions, descriptorHeaders):
        arff_file = self._prepareArff(reactions, descriptorHeaders)
        model_file = self.statsModel.fileName.name

        results_file = "{}_{}.out".format(self.statsModel.pk, uuid.uuid4())
        results_path = os.path.join(settings.TMP_DIR, results_file)

        # Currently, we support only one "response" variable.
        headers = [h for h in reactions.expandedCsvHeaders if h in descriptorHeaders]
        response_index = headers.index(list(self.statsModel.container.outcomeDescriptors)[0].csvHeader) + 1

        # TODO: Validate this input.
        self.wekaPredict(arff_file, model_file, response_index, results_path)

        response = list(self.statsModel.container.outcomeDescriptors)[0]
        if isinstance(response, rxnDescriptors.BoolRxnDescriptor):
            typeConversionFunction = stringToBool
        elif isinstance(response, rxnDescriptors.OrdRxnDescriptor):
            typeConversionFunction = int
        elif isinstance(response, rxnDescriptors.NumRxnDescriptor):
            typeConversionFunction = float
        elif isinstance(response, rxnDescriptors.CatRxnDescriptor):
            typeConversionFunction = str
        else:
            raise TypeError("Response descriptor is of invalid type {}".format(type(respons)))
        results = tuple((reaction, result) for reaction, result in zip(reactions, self._readWekaOutputFile(results_path, typeConversionFunction)))
        return {response: results}


def stringToBool(s):
    if s == 'True':
        return True
    elif s == 'False':
        return False
    else:
        raise ValueError("Tried to convert string to boolean when string was neither 'True' nor 'False' but {}".format(outcome))
