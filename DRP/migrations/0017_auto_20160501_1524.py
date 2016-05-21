# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0016_remove_reaction_legacyid'),
    ]

    operations = [
        migrations.AddField(
            model_name='performedreaction',
            name='legacyID',
            field=models.IntegerField(unique=True, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='performedreaction',
            name='legacyRef',
            field=models.CharField(max_length=40, null=True, blank=True),
        ),
    ]
