# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0005_featureselectioncontainer'),
    ]

    operations = [
        migrations.RenameField(
            model_name='featureselectioncontainer',
            old_name='featureLibrary',
            new_name='featureVisitorLibrary',
        ),
        migrations.RenameField(
            model_name='featureselectioncontainer',
            old_name='featureTool',
            new_name='featureVisitorTool',
        ),
    ]
