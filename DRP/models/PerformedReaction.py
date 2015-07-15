from django.db import models
from Reaction import Reaction
from django.contrib.auth import User

class PerformedReaction(Reaction):
  
  class Meta:
    app_label="DRP"

  user=models.ForeignKey(User)
  performedDateTime=models.DateTimeField('Date Reaction Performed')
  recommendation=models.ForeignKey(RecommendedReaction)
  valid=models.BooleanField()
  public=models.BooleanField()
  duplicateOf=models.ForeignKey(Reaction)

