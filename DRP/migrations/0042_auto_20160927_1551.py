# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import uuid


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0041_auto_20160711_1335'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='boolmoldescriptorvalue',
            name='id',
        ),
        migrations.RemoveField(
            model_name='boolrxndescriptorvalue',
            name='id',
        ),
        migrations.RemoveField(
            model_name='catmoldescriptorvalue',
            name='id',
        ),
        migrations.RemoveField(
            model_name='catrxndescriptorvalue',
            name='id',
        ),
        migrations.RemoveField(
            model_name='nummoldescriptorvalue',
            name='id',
        ),
        migrations.RemoveField(
            model_name='numrxndescriptorvalue',
            name='id',
        ),
        migrations.RemoveField(
            model_name='ordmoldescriptorvalue',
            name='id',
        ),
        migrations.RemoveField(
            model_name='ordrxndescriptorvalue',
            name='id',
        ),
        migrations.AddField(
            model_name='boolmoldescriptorvalue',
            name='uid',
            field=models.CharField(max_length=36, serialize=False, primary_key=True, default=uuid.uuid4),
        ),
        migrations.AddField(
            model_name='boolrxndescriptorvalue',
            name='uid',
            field=models.CharField(max_length=36, serialize=False, primary_key=True, default=uuid.uuid4),
        ),
        migrations.AddField(
            model_name='catmoldescriptorvalue',
            name='uid',
            field=models.CharField(max_length=36, serialize=False, primary_key=True, default=uuid.uuid4),
        ),
        migrations.AddField(
            model_name='catrxndescriptorvalue',
            name='uid',
            field=models.CharField(max_length=36, serialize=False, primary_key=True, default=uuid.uuid4),
        ),
        migrations.AddField(
            model_name='nummoldescriptorvalue',
            name='uid',
            field=models.CharField(max_length=36, serialize=False, primary_key=True, default=uuid.uuid4),
        ),
        migrations.AddField(
            model_name='numrxndescriptorvalue',
            name='uid',
            field=models.CharField(max_length=36, serialize=False, primary_key=True, default=uuid.uuid4),
        ),
        migrations.AddField(
            model_name='ordmoldescriptorvalue',
            name='uid',
            field=models.CharField(max_length=36, serialize=False, primary_key=True, default=uuid.uuid4),
        ),
        migrations.AddField(
            model_name='ordrxndescriptorvalue',
            name='uid',
            field=models.CharField(max_length=36, serialize=False, primary_key=True, default=uuid.uuid4),
        ),
    ]
