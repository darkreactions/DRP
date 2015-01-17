from django.db import models
from django.contrib.auth.models import User

from DRP.data_config import CONFIG

from Lab_Group import Lab_Group
from ModelStats import ModelStats
from Data import Data


class Recommendation(models.Model):
  class Meta:
    app_label = "DRP"

  #Reactant Fields
  for i in CONFIG.reactant_range():
    exec("reactant_{0} = models.CharField(\"Reactant {0}\", max_length=30)".format(i))
    exec("quantity_{0} = models.CharField(\"Quantity {0}\", max_length=10)".format(i))
    exec("unit_{0} = models.CharField(\"Unit {0}\", max_length=4)".format(i))

  score = models.FloatField("Score")
  temp = models.CharField("Temperature", max_length=10)
  time = models.CharField("Time", max_length=10) ###
  pH = models.CharField("pH", max_length=5)

  #Yes/No/? Fields:
  slow_cool = models.CharField("Slow Cool", max_length=10)
  leak = models.CharField("Leak", max_length=10)
  outcome = models.CharField("Outcome", max_length=1)
  purity = models.CharField("Purity", max_length=1)

  #Self-assigning Fields:
  atoms = models.CharField("Atoms", max_length=30, blank=True)
  lab_group = models.ForeignKey(Lab_Group, unique=False)
  model_version = models.ForeignKey(ModelStats, unique=False)
  user = models.ForeignKey(User, unique=False, null=True, blank=True, default=None, related_name="last_user")
  assigned_user = models.ForeignKey(User, unique=False, null=True, blank=True, default=None, related_name="assigned_user")
  seed = models.ForeignKey(Data, unique=False, null=True, blank=True, default=None)
  date_dt = models.DateTimeField("Created", null=True, blank=True)
  complete = models.BooleanField("Complete", default=False)
  seeded = models.BooleanField("From Seed", default=False)

  #Fields for user feedback.
  saved = models.BooleanField("Saved", default=False)
  nonsense = models.BooleanField("Nonsense", default=False)
  hidden = models.BooleanField("Hidden", default=False)
  notes = models.CharField("Notes", max_length=200, blank=True)

  def __unicode__(self):
    return u"REC: {} -- (LAB: {} -- Saved: {})".format(self.score,
                                                       self.lab_group.lab_title,
                                                       self.saved)

