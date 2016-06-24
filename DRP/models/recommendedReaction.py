"""A module containing only the RecommendedReaction class."""
from django.db import models
from .reaction import Reaction


class RecommendedReaction(Reaction):

    """A class to store information about reactions recommended by the machine learning infrastructure."""

    class Meta:
        app_label = 'DRP'

    score = models.FloatField()
    seed = models.ForeignKey(Reaction, null=True, related_name='seeded')
    nonsense = models.BooleanField(default=None)
    hidden = models.BooleanField(default=None)
    saved = models.BooleanField(default=None)
    reference = models.CharField('Text Reference', max_length=200)
