# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0009_auto_20160329_1918'),
    ]

    operations = [
        migrations.RenameField(
            model_name='statsmodel',
            old_name='fileName',
            new_name='outputFile',
        ),
    ]
