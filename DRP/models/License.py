'''A module containing only the License class'''
from django.db import models

class License(models.Model):
  '''A class solely to store the text of a license agreement. These must be signed for someone to use the project site.'''

  class Meta:
    app_label = "DRP"

  text=TextField()
