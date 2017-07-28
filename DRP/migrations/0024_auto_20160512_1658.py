# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0023_auto_20160510_1417'),
    ]

    operations = [
        migrations.AlterField(
            model_name='boolmoldescriptorvalue',
            name='value',
            field=models.NullBooleanField(verbose_name='Value'),
        ),
        migrations.AlterField(
            model_name='boolrxndescriptorvalue',
            name='value',
            field=models.NullBooleanField(verbose_name='Value'),
        ),
    ]
