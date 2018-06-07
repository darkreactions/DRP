# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='ordmoldescriptorvalue',
            name='rater',
        ),
        migrations.RemoveField(
            model_name='ordrxndescriptorvalue',
            name='rater',
        ),
    ]
