from abc import ABCMeta, abstractmethod
import os
from DRP.models import StatsModel, PerformedReaction, TrainingSet, TestSet, TestSetRelation, Descriptor
from django.db.models.fields import AutoField, related
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

      #Copy any non-Foreign Key fields from the descriptor.
      for field in descriptor._meta.fields:
        if not (isinstance(field, AutoField) or
                isinstance(field, related.OneToOneField) or
                isinstance(field, related.ManyToManyField)):
          setattr(pred_descriptor, field.name, getattr(descriptor, field.name) )

      # Add the model's suffix (ID) to the descriptor and name for uniqueness.
      pred_descriptor.heading = descriptor.heading + self._getModelSuffix()
      pred_descriptor.name = descriptor.name + self._getModelSuffix()
      pred_descriptor.model = self.stats_model

      pred_descriptor.save()
      pred_descriptors.append(pred_descriptor)

    self.stats_model.predictsDescriptors.add(*pred_descriptors)
    self.stats_model.save()

  def setSplitter(self, splitter):
    self.stats_model.splitter = splitter.__class__.__name__

  def setModelFile(self, filename):
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
