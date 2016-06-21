# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0004_auto_20160301_0004'),
    ]

    operations = [
        migrations.CreateModel(
            name='FeatureSelectionContainer',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('description', models.TextField(default=b'', blank=True)),
                ('featureLibrary', models.CharField(
                    default=b'', max_length=200, blank=True)),
                ('featureTool', models.CharField(
                    default=b'', max_length=200, blank=True)),
                ('startTime', models.DateTimeField(
                    default=None, null=True, blank=True)),
                ('endTime', models.DateTimeField(
                    default=None, null=True, blank=True)),
                ('built', models.BooleanField(default=False,
                                              verbose_name=b'Has the build procedure been called with this container?', editable=False)),
                ('boolRxnDescriptors', models.ManyToManyField(
                    to='DRP.BoolRxnDescriptor')),
                ('catRxnDescriptors', models.ManyToManyField(
                    to='DRP.CatRxnDescriptor')),
                ('chosenBoolRxnDescriptors', models.ManyToManyField(
                    related_name='chosenForFeatureSelection', to='DRP.BoolRxnDescriptor')),
                ('chosenCatRxnDescriptors', models.ManyToManyField(
                    related_name='chosenForFeatureSelection', to='DRP.CatRxnDescriptor')),
                ('chosenNumRxnDescriptors', models.ManyToManyField(
                    related_name='chosenForFeatureSelection', to='DRP.NumRxnDescriptor')),
                ('chosenOrdRxnDescriptors', models.ManyToManyField(
                    related_name='chosenForFeatureSelection', to='DRP.OrdRxnDescriptor')),
                ('numRxnDescriptors', models.ManyToManyField(
                    to='DRP.NumRxnDescriptor')),
                ('ordRxnDescriptors', models.ManyToManyField(
                    to='DRP.OrdRxnDescriptor')),
                ('outcomeBoolRxnDescriptors', models.ManyToManyField(
                    related_name='outcomeForFeatureSelection', to='DRP.BoolRxnDescriptor')),
                ('outcomeCatRxnDescriptors', models.ManyToManyField(
                    related_name='outcomeForFeatureSelections', to='DRP.CatRxnDescriptor')),
                ('outcomeNumRxnDescriptors', models.ManyToManyField(
                    related_name='outcomeForFeatureSelections', to='DRP.NumRxnDescriptor')),
                ('outcomeOrdRxnDescriptors', models.ManyToManyField(
                    related_name='outcomeForFeatureSelections', to='DRP.OrdRxnDescriptor')),
                ('trainingSet', models.ForeignKey(
                    related_name='trainingSetForFeatureSelection', to='DRP.DataSet', null=True)),
            ],
        ),
    ]
