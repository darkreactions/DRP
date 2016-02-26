# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='recommendedreaction',
            name='hidden',
            field=models.BooleanField(default=None),
        ),
        migrations.AlterField(
            model_name='recommendedreaction',
            name='nonsense',
            field=models.BooleanField(default=None),
        ),
        migrations.AlterField(
            model_name='recommendedreaction',
            name='saved',
            field=models.BooleanField(default=None),
        ),
    ]
