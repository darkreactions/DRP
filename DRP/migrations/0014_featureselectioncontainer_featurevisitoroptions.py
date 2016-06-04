# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0013_auto_20160402_1403'),
    ]

    operations = [
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='featureVisitorOptions',
            field=models.TextField(default=b'{}', blank=True),
        ),
    ]
