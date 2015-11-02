from abc import ABCMeta, abstractmethod
import os
from DRP.models import StatsModel, PerformedReaction, TrainingSet, TestSet, TestSetRelation
from django.conf import settings


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
    reactions = self.getTestingData(self.getModelTag())
    self.predict(reactions, suffix="test")


  # TODO: Change these to use the "relations" schema...
  def setTrainingData(self, reactions):
    for reaction in reactions:
      training_set = TrainingSet(reaction=reaction, model=self.stats_model)
      training_set.save()

  def setTestingData(self, reactions):
    test_set = TestSet()
    test_set.model = self.stats_model
    test_set.name = self.getModelTag()

    # TODO: Need to set in var?
    relations = [TestSetRelation(reaction=reaction, test_set=test_set) for reaction in reactions]

    test_set.save()


  def getTrainingData(self):
    return PerformedReaction.objects.filter(trainingset__model=self.stats_model)

  def getTestingData(self, name):
    return PerformedReaction.objects.filter(testset__name=name)


  def setSplitter(self, splitter):
    self.stats_model.splitter = splitter.__class__.__name__

  def getModelTag(self):
    model = self.stats_model
    return "{}_{}_{}".format(model.library, model.tool, model.id)

  def getModelFilename(self):
    if self.stats_model.fileName:
      return self.stats_model.fileName
    else:
      return os.path.join(settings.MODEL_DIR, self.getModelTag())

