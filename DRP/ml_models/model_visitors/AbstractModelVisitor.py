from abc import ABCMeta, abstractmethod
import os
from DRP.models import StatsModel, PerformedReaction, TrainingSet, TestSet, TestSetRelation, Descriptor
from DRP.models.predRxnDescriptors import PredBoolRxnDescriptor, PredOrdRxnDescriptor, PredNumRxnDescriptor, PredCatRxnDescriptor
from DRP.models.rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, NumRxnDescriptor, CatRxnDescriptor
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, OrdRxnDescriptorValue, NumRxnDescriptorValue, CatRxnDescriptorValue
from django.conf import settings
from django.core.files import File


class AbstractModelVisitor(object):
  __metaclass__ = ABCMeta

  def __init__(self, library, tool, iterations, modelContainer, stats_model=None):

    # Allow this ModelVisitor to "wrap" a pre-existing model instead of building.
    if stats_model is not None:
      if not stats_model.container == modelContainer:
        raise Exception("ModelContainer does not own this stats_model!")
      self.stats_model = stats_model

    # Verify that the modelContainer has been saved and has a valid ID.
    elif modelContainer.id is None:
        raise Exception("ModelContainer object must be saved before making models!")

    else:
      self.stats_model = StatsModel()
      self.stats_model.library = library
      self.stats_model.tool = tool
      self.stats_model.iterations = iterations

      self.stats_model.container = modelContainer
      self.stats_model.save()

    self.DEBUG = False

  def enableDebug(self):
    """ Enables debugging messages. """
    self.DEBUG = True

  @abstractmethod
  def _train(self):
    """A function meant to be overridden by actual ModelVisitor classes.
       The `_train` method should prepare the machine learning model for
       classification and save that model if necessary."""

  @abstractmethod
  def predict(self, reactions, suffix="predict"):
    """Return a dictionary where the key is the response descriptor being predicted
       and the value is an ordered list of predictions for that response where the
       ith prediction corresponds to the ith reaction.

       EG: {<NumRxnDescriptor> "outcome" }:[1,2,1,1]}"""

  def _test(self):
    """A convenience wrapper that gets the model's default testing data,
       predicts outputs from that data, then stores those predictions."""
    reactions = self.getTestingData()
    predictions = self.predict(reactions, suffix="test")
    self.storePredictions(reactions, predictions)

  def summarize(self):
    """Prints the performance metrics of this model using its default test-set."""
    return "Confusion Matrix: {}".format(self.getConfusionMatrices())

  def storePredictions(self, reactions, predictions_dict):
    """Stores the predicted responses in the database as RxnDescriptorValues.
       Specifically, expects the predicted responses as a `predictions_dict`
       in the format specified by the `predict` method."""

    for response, predictions in predictions_dict.items():
      # Get the predictsDescriptor associated with the `response` outcomeDescriptor.
      descriptor = Descriptor.objects.filter(heading=response.heading + self._getModelPredictionSuffix()).downcast().next()

      if isinstance(descriptor, BoolRxnDescriptor):
        val = BoolRxnDescriptorValue()
      elif isinstance(descriptor, OrdRxnDescriptor):
        val = OrdRxnDescriptorValue()
      elif isinstance(descriptor, CatRxnDescriptor):
        val = CatRxnDescriptorValue()
      elif isinstance(descriptor, NumRxnDescriptor):
        val = NumRxnDescriptorValue()
      else:
        error = "Unknown RxnDescriptorValue for '{}'".format(descriptor)
        raise NotImplementedError(error)

      val.descriptor = descriptor
      val.model = self.stats_model

      for reaction, prediction in zip(reactions, predictions):
        # Duplicate the base descriptorValue, but for a new reaction.
        val.id = None
        val.pk = None
        val.value = prediction
        val.reaction = reaction
        val.save()

  def getConfusionMatrices(self):
    """Returns a dicionary of dictionaries of dictionaries, where the outer keys
       are the outcomeDescriptors, the middle keys are the "correct" or "true"
       values, the innermost keys are the "guessed" values that occurred, and
       the value is the integer number of occurrences of that guess when the
       true descriptor was the middle key.

       IE: {field: {true: {guess:#, guess':#},
                    true': {guess:#, guess':#}}
           }
       Eg: {"outcome":
           {"1": {"1": 10
                  "2": 10
                  "3": 13
                  "4": 0
                 }
           , ...
           }
          } """

    return {desc.prediction_of: desc.getConfusionMatrix(self.stats_model)
            for desc in self.stats_model.predictsDescriptors.all().downcast()}


  def setTrainingData(self, reactions):
    """Creates a training-set relation between each provided reaction
       and the stats_model. Note that `reactions` is assumed to be a
       queryset of Reaction objects."""

    # Don't allow more than one TrainingSet to be applied to a single StatsModel
    if TrainingSet.objects.filter(model=self.stats_model).exists():
      raise Exception("This model already has a TrainingSet!")

    for reaction in reactions:
      training_set = TrainingSet(reaction=reaction, model=self.stats_model)
      training_set.save()

  def setTestingData(self, reactions, name=""):
    """Creates a new TestSet full of the provided `reactions` and binds
       that TestSet to the stats_model.

       Optionally applies a unique `name` to that TestSet (instead
       of the default)."""
    test_set = TestSet()
    test_set.model = self.stats_model
    test_set.name = name if name else self.getModelTag()
    test_set.save()

    for reaction in reactions:
      TestSetRelation(reaction=reaction, test_set=test_set).save()

  def setPredictors(self, descriptors):
    """Sets the predictor variables (aka, descriptors) for this model.
       This method assumes all necessary Descriptor objects exist already."""

    if not descriptors.exists():
      raise Exception("Predictor descriptors cannot be empty!")

    self.stats_model.descriptors.add(*descriptors)
    self.stats_model.save()

  def setResponses(self, descriptors):
    """Sets the response variables (aka, outcomeDescriptors) for this model.
       Expects `headings` to be a QuerySet, where each string is the
       heading of an existing Descriptor object.

       Also creates a unique predictsDescriptor entry for each outcomeDescriptor."""

    if not descriptors.exists():
      raise Exception("Response descriptors cannot be empty!")

    self.stats_model.outcomeDescriptors.add(*descriptors)

    self.stats_model.save()

    pred_descriptors = []
    for descriptor in descriptors.downcast():
      if isinstance(descriptor, BoolRxnDescriptor):
        pred_descriptor = PredBoolRxnDescriptor()
      elif isinstance(descriptor, CatRxnDescriptor):
        pred_descriptor = PredCatRxnDescriptor()
      elif isinstance(descriptor, OrdRxnDescriptor):
        pred_descriptor = PredOrdRxnDescriptor()
        pred_descriptor.maximum = descriptor.maximum
        pred_descriptor.minimum = descriptor.minimum
      elif isinstance(descriptor, NumRxnDescriptor):
        pred_descriptor = PredNumRxnDescriptor()
        pred_descriptor.maximum = descriptor.maximum
        pred_descriptor.minimum = descriptor.minimum
      else:
        error = "Predicted Descriptor for '{}'".format(descriptor)
        raise NotImplementedError(error)

      # Add the model's suffix (ID) to the descriptor and name for uniqueness.
      pred_descriptor.heading = descriptor.heading+self._getModelPredictionSuffix()
      pred_descriptor.name = descriptor.name + self._getModelPredictionSuffix()
      pred_descriptor.model = self.stats_model

      pred_descriptor.prediction_of = descriptor
      pred_descriptor.stats_model = self.stats_model

      pred_descriptor.save()
      pred_descriptors.append(pred_descriptor)

    self.stats_model.predictsDescriptors.add(*pred_descriptors)
    self.stats_model.save()

  def setSplitter(self, splitter):
    """Stores the classname of a Splitter object to the stats_model."""
    self.stats_model.splitter = splitter.__class__.__name__

  def setModelFile(self, filename):
    """Uploads a ML-model file to this stats_model."""

    # If a file is already in the models directory, don't re-upload it.
    destined_path = os.path.join(settings.MODEL_DIR, filename)
    if os.path.isfile(destined_path):
      self.stats_model.fileName.name = filename
      self.stats_model.save()

    else:
      if self.DEBUG:
        print "Uploading model file to: {}".format(destined_path)
      f = open(filename, 'r')
      self.stats_model.fileName.save(filename, File(f))


  def getTrainingData(self):
    """Returns a queryset of the reactions used to train this ML-model."""
    return PerformedReaction.objects.filter(trainingset__model=self.stats_model)

  def getTestingData(self, testset_name=""):
    """Returns a queryset of the reactions used to test this ML-model.

       Optionally accepts the name of a specific testset to retrieve."""
    if not testset_name: testset_name = self.getModelTag()
    return PerformedReaction.objects.filter(testset__name=testset_name)

  def _getModelPredictionSuffix(self):
    """Returns a unique "suffix" for the predictions of this stats_model."""
    return "_predicted_{}".format(self.stats_model.id)

  def getPredictors(self):
    """Returns a queryset of descriptors used by this stats_model for
       producing predictions of the various response variables."""
    return self.stats_model.descriptors.all()

  def getResponses(self):
    """Returns a queryset of outcomeDescriptors used by this stats_model."""
    return self.stats_model.outcomeDescriptors.all()

  def getModelTag(self):
    """Returns a unique "name" for this stats_model."""
    container = self.stats_model.container
    return "{}_{}_{}".format(container.library, container.tool, self.stats_model.id)

  def getModelFilename(self):
    """Returns the filename of the ML-model file."""
    if self.stats_model.fileName.name:
      return self.stats_model.fileName.name
    else:
      raise IOError("No file has been uploaded to this stats_model.")
