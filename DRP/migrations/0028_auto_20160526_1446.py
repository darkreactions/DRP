# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0027_auto_20160526_1443'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compoundquantity',
            name='amount',
            field=models.DecimalField(decimal_places=5, validators=[DRP.models.validators.GreaterThanValidator(
                0)], max_digits=12, blank=True, help_text='(in mmoles, 5 decimal places)', null=True),
        ),
    ]
