# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0011_statsmodel_inputfile'),
    ]

    operations = [
        migrations.AddField(
            model_name='modelcontainer',
            name='modelVisitorOptions',
            field=models.TextField(default=b'', blank=True),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='splitterOptions',
            field=models.TextField(default=b'', blank=True),
        ),
    ]
