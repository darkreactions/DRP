# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0040_auto_20160707_1515'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compoundquantity',
            name='amount',
            field=models.DecimalField(null=True, max_digits=12, decimal_places=5, help_text='(in mmoles, up to 5 decimal places)', blank=True, validators=[DRP.models.validators.GreaterThanValidator(0)]),
        ),
        migrations.AlterField(
            model_name='labgroup',
            name='defaultDescriptors',
            field=models.ManyToManyField(related_name='isDefaultForLabGroups', to='DRP.Descriptor', blank=True),
        ),
    ]
