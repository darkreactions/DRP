"""A module containing only the License clas."""
from django.db import models


class License(models.Model):

    """A class solely to store the text of a license agreement. These must be signed for someone to use the project site."""

    class Meta:
        app_label = "DRP"
        get_latest_by = 'effectiveDate'

    text = models.TextField()
    effectiveDate = models.DateField(verbose_name='Effective Date', help_text='The license will become effective on midnight of the provided date.')
