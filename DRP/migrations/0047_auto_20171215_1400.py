# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('DRP', '0046_auto_20171215_1337'),
    ]

    operations = [
        migrations.AddField(
            model_name='ordmoldescriptorvalue',
            name='rater',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, default=None),
        ),
        migrations.AddField(
            model_name='ordrxndescriptorvalue',
            name='rater',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, default=None),
        ),
    ]
