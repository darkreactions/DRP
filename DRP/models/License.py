from django.db import models

class License(models.Model):

  class Meta:
    app_label = "DRP"

  text=TextField()
