# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0031_auto_20160620_1357'),
    ]

    operations = [
        migrations.AddField(
            model_name='compound',
            name='labGroups',
            field=models.ManyToManyField(to='DRP.LabGroup', verbose_name=b'Lab Groups'),
        ),
    ]
