from django.conf import settings
import uuid
from DRP.ml_models.feature_visitors.AbstractFeatureVisitor import AbstractFeatureVisitor, logger
from django.core.exceptions import ImproperlyConfigured
import subprocess
import os
from abc import abstractmethod


class AbstractWekaFeatureVisitor(AbstractFeatureVisitor):

    maxResponseCount = 1
    
    def __init__(self, container, *args, **kwargs):
        super(AbstractWekaFeatureVisitor, self).__init__(*args, **kwargs)

        self.container = container
        
        self.WEKA_VERSION = "3.6" # The version of WEKA to use.
    
    def _prepareArff(self, reactions, whitelistHeaders, verbose=False):
        """Writes an *.arff file using the provided queryset of reactions."""
        logger.debug("Preparing ARFF file...")
        filename = "featureSelection_{}_{}.arff".format(self.container.pk, uuid.uuid4())
        filepath = os.path.join(settings.TMP_DIR, filename)
        while os.path.isfile(filepath): #uber paranoid making sure we don't race condition
            filename = "{}_{}.arff".format(self.container.pk, uuid.uuid4())
            filepath = os.path.join(settings.TMP_DIR, filename)
        if verbose:
            print "Writing arff to {}".format(filepath)
        with open(filepath, "w") as f:
          reactions.toArff(f, expanded=True, whitelistHeaders=whitelistHeaders)
        return filepath

    def _readWekaOutput(self, output):
        """
        Reads a weka feature selection output and outputs a list of descriptors
        """
        start_line = "Selected attributes:"
        raw_lines = output.split('\n')
        
        found = False
        descriptors = []
        for line in raw_lines:
            if found:
                desc = line.strip()
                if desc: #check if line is just whitespace
                    descriptors.append(desc)
            if line.startswith(start_line):
                found = True

        return descriptors


    def _runWekaCommand(self, command, verbose=False):
        """Sets the CLASSPATH necessary to use Weka, then runs a shell `command`."""
        if not settings.WEKA_PATH[self.WEKA_VERSION]:
            raise ImproperlyConfigured("'WEKA_PATH' is not set in settings.py!")
        set_path = "export CLASSPATH=$CLASSPATH:{}; ".format(settings.WEKA_PATH[self.WEKA_VERSION])
        command = set_path + command
        logger.debug("Running in Shell:\n{}".format(command))
        if verbose:
            print "Running in Shell:\n{}".format(command)
        output = subprocess.check_output(command, shell=True)
        return output
    
    def train(self, reactions, descriptorHeaders, verbose=False):
        arff_file = self._prepareArff(reactions, descriptorHeaders, verbose=verbose)

        results_file = "featureSelection_{}_{}.out".format(self.container.pk, uuid.uuid4())
        results_path = os.path.join(settings.TMP_DIR, results_file)
        
        # Currently, we support only one "response" variable.
        headers = [h for h in reactions.expandedCsvHeaders if h in descriptorHeaders]
        response_index = headers.index(list(self.container.outcomeDescriptors)[0].csvHeader) + 1
        command = self.wekaTrainCommand(arff_file, response_index)
        
        output = self._runWekaCommand(command, verbose=verbose)
        if verbose:
            print output

        descriptor_headers = self._readWekaOutput(output)
        return descriptor_headers
        
        

    @abstractmethod
    def wekaTrainCommand(self, arff_file, response_index, out_file):
        """Abstract method that returns the weka train command"""
