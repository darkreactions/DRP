# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0038_auto_20160630_1317'),
    ]

    operations = [
        migrations.AlterField(
            model_name='reaction',
            name='labGroup',
            field=models.ForeignKey(verbose_name='Lab Group', to='DRP.LabGroup'),
        ),
    ]
