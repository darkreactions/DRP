# -*- coding: utf-8 -*-
# You may need to fake this migration
# sorry :-(
from __future__ import unicode_literals

from django.db import migrations, models
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0007_delete_compound_quantity_manager'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='modelcontainer',
            name='featureLibrary',
        ),
    ]
