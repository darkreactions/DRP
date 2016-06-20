# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0030_auto_20160605_2127'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compound',
            name='labGroup',
            field=models.ForeignKey(related_name='old_compounds', verbose_name=b'Lab Group', to='DRP.LabGroup'),
        ),
    ]
