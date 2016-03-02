from django.conf import settings
import uuid
from DRP.ml_models.feature_visitors.AbstractFeatureVisitor import AbstractFeatureVisitor, logger
from django.core.exceptions import ImproperlyConfigured
import subprocess
import os
from abc import abstractmethod


class AbstractWekaFeatureVisitor(AbstractFeatureVisitor):

    maxResponseCount = 1
    
    def __init__(self, *args, **kwargs):
    super(AbstractWekaFeatureVisitor, self).__init__(*args, **kwargs)
    
    self.WEKA_VERSION = "3.6" # The version of WEKA to use.
    
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
    
