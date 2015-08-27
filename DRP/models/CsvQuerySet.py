'''A module for turning querysets into csv files or output'''
import csv
from django.db import models
import abc

class CsvModel(models.Model):

  class Meta:
    abstract = True

  @abc.abstractproperty
  def values(self):
    '''Returns a dict of values suitable for use by a csv.DictWriter'''
    return {field.name:getattr(self, field.name) for field in self._meta.fields}

  @abc.abstractproperty
  def expandedValues(self):
    '''Returns a dict of values suitable for use by a csv.DictWriter, including any special expanded values'''
    return self.values

class CsvQuerySet(models.query.QuerySet):
  '''This queryset permits the output of the data from a model as a csv'''

  __metaclass__ = abc.ABCMeta

  @abc.abstractproperty
  def expandedHeaders(self):
    '''For some classes, like Compounds and Reactions, there needs to be the option to send additional data (descriptors)
    as part of the csv. This method permits that expansion, and defaults to return the non-expanded headers.'''
    return self.headers

  @abc.abstractproperty
  def headers(self):
    '''The basic headers to be used for the model. Note that the implementation on the CsvQuerySet class is extremely basic,
    and will fail if any field holds a relationship, and will not include automagically generated fields.'''
    return [field.name for field in self.model._meta.fields]

  def toCsv(self, writeable, expanded=False, missing="?"):#TODO:figure out most sensible default for missing values
    '''Writes the csv data to the writeable (file, or for Django a HttpResponse) object. Expanded outputs any expanded
      information that the corresponding methods provide- this requires the model being called to have a property
      'expandedValues', which should be a dictionary like object of values, using fieldNames as keys as output
      by fetchExpandedHeaders.
    '''

    if expanded:
      writer = csv.DictWriter(writeable, fieldnames=self.expandedHeaders, restval=missing)
    else:
      writer = csv.DictWriter(writeable, fieldnames=self.headers, restval=missing)

    writer.writeheader()
    for item in self:
      if expanded:
        writer.writerow(item.expandedValues)
      else:
        writer.writerow(item.values)
