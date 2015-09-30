from django.db import models

class CompoundRole(models.Model):

  class Meta:
    verbose_name='Compound Role Category'
    verbose_name_plural='Compound Role Categories'
    app_label='DRP'

  label=models.CharField(max_length=255, unique=True)
  description=models.TextField()

  def __unicode__(self):
    return self.label
