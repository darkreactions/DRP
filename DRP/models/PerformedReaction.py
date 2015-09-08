'''A module containing only the PerformedReaction class'''
from django.db import models
from Reaction import Reaction
from RecommendedReaction import RecommendedReaction
from django.contrib.auth.models import User
from StatsModel import StatsModel

class PerformedReaction(Reaction):
  '''A class representing concrete instances of reactions that have actually been performed'''
  
  class Meta:
    app_label="DRP"

  user=models.ForeignKey(User)
  reference=models.CharField(unique=True, max_length=40) 
  performedDateTime=models.DateTimeField('Date Reaction Performed')
  recommendation=models.ForeignKey(RecommendedReaction, unique=False, null=True, default=None, related_name='resultantExperiment')
  legacyRecommendedFlag=models.NullBooleanField(default=None)
  '''If this reaction was based from a recommendation, reference that recommendation'''
  valid=models.BooleanField()
  '''A flag to denote reactions which have been found to be invalid, for instance,
  if the wrong reactant was used or some bad lab record has been found'''
  public=models.BooleanField()
  '''We are storing this for legacy purposes but not using it in this version'''
  duplicateOf=models.ForeignKey('self', related_name='duplicatedBy')
  '''Describes the many to many mapping when a StatsModel uses a Performed Reaction as part of its
  test or training sets. Must be placed on this model as a workaround to circular dependency issues'''
  inTrainingSetFor=models.ManyToManyField(StatsModel, related_name='trainingSet')
  inTestSetFor=models.ManyToManyField(StatsModel, related_name='testSet')
