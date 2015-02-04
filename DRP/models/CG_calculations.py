from django.db import models

class CG_calculations(models.Model):
  class Meta:
    app_label = "DRP"

  json_data = models.TextField()
  compound = models.CharField(max_length=200)
  smiles = models.CharField(max_length=255, unique=True)
  json = models.TextField(null=True, default="{}")

  def __unicode__(self):
    return u"{} ({})".format(self.compound, self.smiles)

