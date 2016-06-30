# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0003_auto_20160229_2044'),
    ]

    operations = [
        migrations.AlterField(
            model_name='metriccontainer',
            name='description',
            field=models.TextField(default='', blank=True),
        ),
        migrations.AlterField(
            model_name='metriccontainer',
            name='endTime',
            field=models.DateTimeField(default=None, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='metriccontainer',
            name='startTime',
            field=models.DateTimeField(default=None, null=True, blank=True),
        ),
    ]
