from django.db import models

class Descriptor(models.Model):
  
  class Meta:
    app_label='DRP'

  heading=models.CharField(max_length=200)
  name=models.CharField('Full name', max_length=300)
  kind=models.Charfield('Kind', max_length=20, choices=(('Cat', 'Categorical'), ('Ord', 'Ordinal'), ('Num', 'Numerical'), ('Bool', 'Boolean'))
  NaNCat = models.CharField('Not-a-value placeholder for the case of a categorical descriptor', max_length=200, null=True, default=None)
  NaNOrd = models.CharField('Not-a-value placeholder for the case of an ordinal descriptor', max_length=100, null=True, default=None)
  NaNNum = models.CharField('Not-a-value placeholder for the case of a numerical descriptor', max_length=100, null=True, default=None)
  NaNBool = models.CharField('Not-a-value placeholder for the case of a boolean descriptor', max_length=100, null=True, default=None)
