# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0014_featureselectioncontainer_featurevisitoroptions'),
    ]

    operations = [
        migrations.AddField(
            model_name='reaction',
            name='legacyID',
            field=models.IntegerField(unique=True, null=True, blank=True),
        ),
    ]
