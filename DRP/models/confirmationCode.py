"""Confirmation code module used in email validation."""
from django.db import models
from django.contrib.auth.models import User


class ConfirmationCode(models.Model):
    """
    Class to store confirmation codes for emails.

    This class stores confirmation codes for users so that we can confirm that they have control over
    Their provided email address- this ensures that the license agreement they sign is valid.
    """

    class Meta:
        app_label = "DRP"

    user = models.OneToOneField(User)
    code = models.CharField(max_length=36, unique=True)
