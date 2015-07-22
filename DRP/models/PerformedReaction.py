'''A module containing only the PerformedReaction class'''
from django.db import models
from Reaction import Reaction
from django.contrib.auth import User

class PerformedReaction(Reaction):
  '''A class representing concrete instances of reactions that have actually been performed'''
  
  class Meta:
    app_label="DRP"

  user=models.ForeignKey(User)
  performedDateTime=models.DateTimeField('Date Reaction Performed')
  recommendation=models.ForeignKey(RecommendedReaction)
  '''If this reaction was based from a recommendation, reference that recommendation'''
  valid=models.BooleanField()
  '''A flag to denote reactions which have been found to be invalid, for instance,
  if the wrong reactant was used or some bad lab record has been found'''
  public=models.BooleanField()
  duplicateOf=models.ForeignKey(Reaction)

