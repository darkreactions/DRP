from django.db import models

class StatsModelTags(models.Model)

  class Meta:
    app_label='DRP'

  text=models.Charfield(max_length=200, unique=True)
