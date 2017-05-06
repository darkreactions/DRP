# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0043_auto_20160927_1553'),
    ]

    operations = [
        migrations.AddField(
            model_name='compound',
            name='calculating',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='compound',
            name='dirty',
            field=models.BooleanField(default=True),
        ),
        migrations.AddField(
            model_name='compound',
            name='recalculate',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='reaction',
            name='calculating',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='reaction',
            name='dirty',
            field=models.BooleanField(default=True),
        ),
        migrations.AddField(
            model_name='reaction',
            name='recalculate',
            field=models.BooleanField(default=False),
        ),
    ]
