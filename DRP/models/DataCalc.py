from django.db import models
import json

class DataCalc(models.Model):
  class Meta:
    app_label = "DRP"

  contents = models.TextField(default="{}")

  def __unicode__(self):
    return u"{}".format(self.contents);

  #Convert the stringy contents to an actual array/JSON object.
  def make_json(self):
    return json.loads(self.contents)

