"""A compound role signifies the role a compound plays in a specific reaction."""
from django.db import models


class CompoundRole(models.Model):
    """
    A compound role signifies the role a compound plays in a specific reaction.

    This is included as a part of the ontology by being a component of
    a compound quantity.
    """

    class Meta:
        verbose_name = 'Compound Role Category'
        verbose_name_plural = 'Compound Role Categories'
        app_label = 'DRP'

    label = models.CharField(max_length=255, unique=True)
    description = models.TextField()

    def __str__(self):
        """Return the label string for the unicode rep."""
        return self.label
