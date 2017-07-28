# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.db import models, migrations
import uuid
import time

runset = set()
def myuid(apps, schema_editor):
    Models = [apps.get_model('DRP', 'BoolRxnDescriptorValue'), 
            apps.get_model('DRP', 'BoolMolDescriptorValue'),
            apps.get_model('DRP', 'CatMolDescriptorValue'),
            apps.get_model('DRP', 'CatRxnDescriptorValue'),
            apps.get_model('DRP', 'NumRxnDescriptorValue'),
            apps.get_model('DRP', 'NumMolDescriptorValue'),
            apps.get_model('DRP', 'OrdRxnDescriptorValue'),
            apps.get_model('DRP', 'OrdMolDescriptorValue'),
    ]
    for m in Models:
        for o in m.objects.all():
            o.uid=uuid.uuid4() #this probably could use some refinement...
            o.save()


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0041_auto_20160711_1335'),
    ]

    operations = [
        migrations.AddField(
            model_name='boolmoldescriptorvalue',
            name='uid',
            field=models.CharField(
                max_length=36, serialize=False, null=True),
        ),
        migrations.AddField(
            model_name='boolrxndescriptorvalue',
            name='uid',
            field=models.CharField(
                max_length=36, serialize=False, null=True),
        ),
        migrations.AddField(
            model_name='catmoldescriptorvalue',
            name='uid',
            field=models.CharField(
                max_length=36, serialize=False, null=True),
        ),
        migrations.AddField(
            model_name='catrxndescriptorvalue',
            name='uid',
            field=models.CharField(
                max_length=36, serialize=False, null=True),
        ),
        migrations.AddField(
            model_name='nummoldescriptorvalue',
            name='uid',
            field=models.CharField(
                max_length=36, serialize=False, null=True),
        ),
        migrations.AddField(
            model_name='numrxndescriptorvalue',
            name='uid',
            field=models.CharField(
                max_length=36, serialize=False, null=True),
        ),
        migrations.AddField(
            model_name='ordmoldescriptorvalue',
            name='uid',
            field=models.CharField(
                max_length=36, serialize=False, null=True),
        ),
        migrations.AddField(
            model_name='ordrxndescriptorvalue',
            name='uid',
            field=models.CharField(
                max_length=36, serialize=False, null=True),
        ),
        migrations.RunPython(myuid, reverse_code=migrations.RunPython.noop)
    ]
