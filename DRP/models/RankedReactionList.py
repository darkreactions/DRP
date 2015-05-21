from django.db import models

from django.contrib.auth.models import User

import json

#TODO: Are ranked reactions even necessary? Probably not... -- Casey (1/16/15)
class RankedReactionList(models.Model):
  class Meta:
    app_label = "DRP"

  original_list = models.TextField()
  seed = models.TextField()
  ranked_list = models.TextField()
  ranker = models.ForeignKey(User, unique=False, default=None, null=True)

  def __unicode__(self):
    return u"RANKED_LIST: Seed: {} -- (Ranker: {})".format(self.seed, self.ranker)

  def get_original_list(self):
    return json.loads(self.original_list)

  def get_ranked_list(self):
    return json.loads(self.ranked_list)

  def get_seed(self):
    return json.loads(self.seed)

  #Returns a shuffled list and the index of the seed.
  def get_shuffled_list(self):
    shuffled_list = self.get_original_list()
    random.shuffle(shuffled_list)
    return (shuffled_list)

  def store_ranked_list(self, ranked_list, ranker):
    self.ranker = ranker
    self.ranked_list = json.dumps(ranked_list)
    self.save()



