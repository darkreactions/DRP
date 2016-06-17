# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0021_auto_20160501_1815'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compoundquantity',
            name='amount',
            field=models.FloatField(
                help_text=b'(in mmoles)', null=True, blank=True),
        ),
    ]
