from abc import ABCMeta, abstractmethod
import os, time
from DRP.models import StatsModel, PerformedReaction, TrainingSet, TestSet, TestSetRelation, Descriptor
from DRP.models.rxnDescriptorValues import getDescriptorAndEmptyVal
from django.conf import settings
from django.core.files import File


class AbstractModelVisitor(object):
  __metaclass__ = ABCMeta

  def __init__(self):
    self.stats_model = StatsModel()
    self.DEBUG = False

  def enableDebug(self):
    self.DEBUG = True

  @abstractmethod
  def _train(self):
    pass

  @abstractmethod
  def predict(self, reactions, suffix="predict"):
    pass

  def _test(self):
    reactions = self.getTestingData()
    predictions = self.predict(reactions, suffix="test")
    self.storePredictions(reactions, predictions)

  def storePredictions(self, reactions, predictions):
    for response in self.getResponses():
      desc, val = getDescriptorAndEmptyVal(response.heading + self._getModelSuffix())
      val.descriptor = desc
      val.model = self.stats_model

      for reaction, prediction in zip(reactions, predictions):
        # Duplicate the base descriptorValue, but for a new reaction.
        val.id = None
        val.pk = None
        val.value = prediction
        val.reaction = reaction
        val.save()


  def setTrainingData(self, reactions):
    for reaction in reactions:
      training_set = TrainingSet(reaction=reaction, model=self.stats_model)
      training_set.save()

  def setTestingData(self, reactions):
    test_set = TestSet()
    test_set.model = self.stats_model
    test_set.name = self.getModelTag()
    test_set.save()

    for reaction in reactions:
      TestSetRelation(reaction=reaction, test_set=test_set).save()

  def setPredictors(self, headers):
    descriptors = [Descriptor.objects.get(heading=header) for header in headers]
    self.stats_model.descriptors.add(*descriptors)
    self.stats_model.save()

  def setResponses(self, headers):
    descriptors = [Descriptor.objects.get(heading=h).downcast() for h in headers]
    self.stats_model.outcomeDescriptors.add(*descriptors)

    self.stats_model.save()

    pred_descriptors = []
    for descriptor in descriptors:
      # Copy the descriptor to a pred_descriptor so we retain the descriptor type.
      pred_descriptor = descriptor.__class__()
      # TODO: Copy all fields from descriptor into pred_descriptor
      pred_descriptor.heading = descriptor.heading + self._getModelSuffix()
      pred_descriptor.name = descriptor.name + self._getModelSuffix()
      pred_descriptor.maximum = descriptor.maximum #TODO: HARD CODED
      pred_descriptor.minimum = descriptor.minimum #TODO: HARD CODED
      pred_descriptor.model = self.stats_model
      pred_descriptor.save()
      pred_descriptors.append(pred_descriptor)

    self.stats_model.predictsDescriptors.add(*pred_descriptors)
    self.stats_model.save()

  def setSplitter(self, splitter):
    self.stats_model.splitter = splitter.__class__.__name__

  def setModelFile(self, filename):
    f = open(filename, 'r')
    self.stats_model.fileName.save(filename, File(f))

    # Wait for the file to upload.
    max_wait = 10 # seconds
    filepath = os.path.join(settings.BASE_DIR, settings.MODEL_DIR, filename)

    while max_wait > 0:
      try:
        with open(filepath, 'rb'):
          break
      except IOError:
        if self.DEBUG:
          print "Waiting {} for file to be written: {}".format(max_wait, filename)
        time.sleep(1)
        max_wait -= 1


  def getTrainingData(self):
    return PerformedReaction.objects.filter(trainingset__model=self.stats_model)

  def getTestingData(self):
    return PerformedReaction.objects.filter(testset__name=self.getModelTag())

  def _getModelSuffix(self):
    return "_{}".format(self.stats_model.id)

  def getPredictors(self):
    return self.stats_model.descriptors.all()

  def getResponses(self):
    return self.stats_model.outcomeDescriptors.all()

  def getModelTag(self):
    model = self.stats_model
    return "{}_{}_{}".format(model.library, model.tool, model.id)

  def getModelFilename(self):
    return self.stats_model.fileName.name
