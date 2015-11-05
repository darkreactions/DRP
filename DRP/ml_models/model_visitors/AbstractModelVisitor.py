from abc import ABCMeta, abstractmethod
import os
from DRP.models import StatsModel, PerformedReaction, TrainingSet, TestSet, TestSetRelation, Descriptor
from django.conf import settings
from django.core.files import File


class AbstractModelVisitor(object):
  __metaclass__ = ABCMeta

  def __init__(self):
    self.stats_model = StatsModel()
    self.debug = False

  def enableDebug(self):
    self.debug = True

  @abstractmethod
  def _train(self):
    pass

  @abstractmethod
  def predict(self, reactions, suffix="predict"):
    pass

  def _test(self):
    reactions = self.getTestingData()
    self.predict(reactions, suffix="test")

  # TODO: Store predictions as necessary.
  def storePredictions(reactions, predictions, responseField):
    pass

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


  def getTrainingData(self):
    return PerformedReaction.objects.filter(trainingset__model=self.stats_model)

  def getTestingData(self):
    return PerformedReaction.objects.filter(testset__name=self.getModelTag())

  def setPredictors(self, headers):
    descriptors = [Descriptor.objects.get(heading=header) for header in headers]
    self.stats_model.descriptors.add(*descriptors)
    self.stats_model.save()

  def setResponses(self, headers):
    descriptors = [Descriptor.objects.get(heading=header) for header in headers]
    self.stats_model.outcomeDescriptors.add(*descriptors)
    self.stats_model.save()

  def getPredictors(self):
    return self.stats_model.descriptors.all()

  def getResponses(self):
    return self.stats_model.outcomeDescriptors.all()

  def setSplitter(self, splitter):
    self.stats_model.splitter = splitter.__class__.__name__

  def getModelTag(self):
    model = self.stats_model
    return "{}_{}_{}".format(model.library, model.tool, model.id)

  def getModelFilename(self):
    if self.stats_model.fileName:
      return self.stats_model.fileName.name
    else:
      filename = "{}.model".format(self.getModelTag())
      return os.path.join(settings.MODEL_DIR, filename)

  def setModelFile(self, filename):
    f = open(filename, 'r')
    self.stats_model.fileName.save(filename, File(f))

