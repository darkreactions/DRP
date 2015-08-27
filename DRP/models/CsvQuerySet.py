'''A module for turning querysets into csv files or output'''
import csv
from django.db import models

class CsvQuerySet(models.query.QuerySet):
  '''This queryset permits the output of the data from a model as a csv'''

  @property
  def expandedHeaders(self):
    '''For some classes, like Compounds and Reactions, there needs to be the option to send additional data (descriptors)
    as part of the csv. This method permits that expansion, and defaults to have no effect.'''
    return []

  def toCsv(self, writeable, expanded=False, NaN='NaN')#TODO:figure out most sensible default
    '''Writes the csv data to the writeable (file, or for Django a HttpResponse) object. Expanded outputs any expanded
      information that the corresponding methods provide- this requires the model being called to have a property
      'expandedValues', which should be a dictionary like object of values, using fieldNames as keys as output
      by fetchExpandedHeaders.
    '''

    baseFieldNames = self.model_meta.get_all_field_names()
    exFieldNames = self.expandedHeaders

    writer = csv.writer(writeable):
    writer.writerow(baseFieldNames+exFieldNames)
    for item in self:
      if expanded:
        writer.writerow([NaN if getattr(item, fieldName) is None else getattr(item, fieldName) for fieldName in baseFieldNames] + [NaN if item.expandedValues.get(fieldName) is None else item.expandedValues.get(fieldName) for fieldName in exFieldNames])
