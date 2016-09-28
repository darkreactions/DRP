# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.molDescriptorValues
import DRP.models.rxnDescriptorValues


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0042_auto_20160927_1551'),
    ]

    operations = [
        migrations.AlterField(
            model_name='boolmoldescriptorvalue',
            name='uid',
            field=models.CharField(default=DRP.models.molDescriptorValues.molUid, max_length=36, serialize=False, primary_key=True),
        ),
        migrations.AlterField(
            model_name='boolrxndescriptorvalue',
            name='uid',
            field=models.CharField(default=DRP.models.rxnDescriptorValues.rxnUid, max_length=36, serialize=False, primary_key=True),
        ),
        migrations.AlterField(
            model_name='catmoldescriptorvalue',
            name='uid',
            field=models.CharField(default=DRP.models.molDescriptorValues.molUid, max_length=36, serialize=False, primary_key=True),
        ),
        migrations.AlterField(
            model_name='catrxndescriptorvalue',
            name='uid',
            field=models.CharField(default=DRP.models.rxnDescriptorValues.rxnUid, max_length=36, serialize=False, primary_key=True),
        ),
        migrations.AlterField(
            model_name='nummoldescriptorvalue',
            name='uid',
            field=models.CharField(default=DRP.models.molDescriptorValues.molUid, max_length=36, serialize=False, primary_key=True),
        ),
        migrations.AlterField(
            model_name='numrxndescriptorvalue',
            name='uid',
            field=models.CharField(default=DRP.models.rxnDescriptorValues.rxnUid, max_length=36, serialize=False, primary_key=True),
        ),
        migrations.AlterField(
            model_name='ordmoldescriptorvalue',
            name='uid',
            field=models.CharField(default=DRP.models.molDescriptorValues.molUid, max_length=36, serialize=False, primary_key=True),
        ),
        migrations.AlterField(
            model_name='ordrxndescriptorvalue',
            name='uid',
            field=models.CharField(default=DRP.models.rxnDescriptorValues.rxnUid, max_length=36, serialize=False, primary_key=True),
        ),
    ]
