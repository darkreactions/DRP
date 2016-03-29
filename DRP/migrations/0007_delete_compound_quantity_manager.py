# -*- coding: utf-8 -*-
# Something really weird happened with the CompoundQuantityManager,
# which was supposed to be a manager but accidentally inherited from model
# I think it has to do with the transition from South
# You may need to fake this migration to proceed

from __future__ import unicode_literals

from django.db import migrations, models
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0006_auto_20160302_2159'),
    ]

    operations = [
        migrations.DeleteModel(
            name='CompoundQuantityManager',
        ),
    ]
