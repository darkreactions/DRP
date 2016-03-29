# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0007_validateCalcSoftware_removeFeatureTools'),
    ]

    operations = [
        migrations.DeleteModel(
           name='CompoundQuantityManager',
        ),
    ]
