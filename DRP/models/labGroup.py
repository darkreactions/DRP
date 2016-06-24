"""A module containing only the LabGroup class."""
from django.db import models
from django.contrib import auth
from django.contrib.auth.models import User
from django.contrib.auth.hashers import make_password
from django.conf import settings


class LabGroupManager(models.Manager):

    """A custom manager with a convenience function so that we can create new lab groups."""

    def makeLabGroup(self, title, address, email, access_code):
        """A function to create new lab groups easily."""
        return LabGroup(title=title, address=address, email=email, access_code=make_password(access_code, settings.LAB_GROUP_HASH_SALT))


class LabGroup(models.Model):

    """A class for describing a collection of scientists belonging to the same group."""

    class Meta:
        app_label = "DRP"
        verbose_name = 'Lab Group'

    title = models.CharField(max_length=200, unique=True, error_messages={
                             'unique': "This name is already taken."})
    address = models.CharField(max_length=200)
    email = models.CharField(max_length=254,  # Maximum length of email address
                             default='')
    access_code = models.CharField(max_length=128)
    """An access code to allow members to \'prove\' that they are permitted to join the lab. Normally held by a laboratory administrator."""
    legacy_access_code = models.CharField(max_length=20)
    """An older version of the access code. Made a part of this model for legacy support."""
    users = models.ManyToManyField(User, blank=True)
    objects = LabGroupManager()

    def __unicode__(self):
        """Use the title as the unicode rep."""
        return self.title
