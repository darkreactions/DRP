# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0026_auto_20160526_1429'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compoundquantity',
            name='amount',
            field=models.FloatField(blank=True, help_text=b'(in mmoles)', null=True, validators=[
                                    DRP.models.validators.GreaterThanValidator(0)]),
        ),
    ]
