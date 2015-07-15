from django.db import models

class DescriptorValue:

  class Meta:
    app_label="DRP"

  descriptor = models.ForeignKey(Descriptor)
  booleanValue= models.BooleanField('Value if descriptor is a boolean', null=True)
  ordValue = models.PositiveIntegerField('Value if descriptor is an ordinal', null=True)
  catValue = models.Charfield('Value if descriptor is a category', max_length=200, null=True)
  numValue = models.FloatField('Value if descritpor is continuous', null=True)
  isPredicted=models.BooleanField()
