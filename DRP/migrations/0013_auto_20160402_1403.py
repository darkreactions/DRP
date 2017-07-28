# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0012_auto_20160402_1319'),
    ]

    operations = [
        migrations.AlterField(
            model_name='modelcontainer',
            name='modelVisitorOptions',
            field=models.TextField(default='{}', blank=True),
        ),
        migrations.AlterField(
            model_name='modelcontainer',
            name='splitterOptions',
            field=models.TextField(default='{}', blank=True),
        ),
    ]
