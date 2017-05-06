# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0035_auto_20160621_1057'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='compound',
            name='abbrev',
        ),
        migrations.RemoveField(
            model_name='compound',
            name='labGroup',
        ),
        migrations.AlterField(
            model_name='compound',
            name='CSID',
            field=models.PositiveIntegerField(
                unique=True, null=True, verbose_name='Chemspider ID'),
        ),
    ]
