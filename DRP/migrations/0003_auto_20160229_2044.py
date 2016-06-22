# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='metriccontainer',
            name='built',
            field=models.BooleanField(
                default=False, verbose_name=b'Has the build procedure been called with this container?', editable=False),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='description',
            field=models.TextField(default=b''),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='endTime',
            field=models.DateTimeField(default=None, null=True),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='fileName',
            field=models.FileField(
                max_length=200, upload_to=b'metrics', blank=True),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='invalid',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='startTime',
            field=models.DateTimeField(default=None, null=True),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='trainingSet',
            field=models.ForeignKey(
                related_name='trainingSetForMetric', to='DRP.DataSet', null=True),
        ),
        migrations.AlterField(
            model_name='modelcontainer',
            name='description',
            field=models.TextField(blank=True),
        ),
        migrations.AlterField(
            model_name='modelcontainer',
            name='featureLibrary',
            field=models.CharField(default=b'', max_length=200, blank=True),
        ),
        migrations.AlterField(
            model_name='modelcontainer',
            name='featureTool',
            field=models.CharField(default=b'', max_length=200, blank=True),
        ),
        migrations.AlterField(
            model_name='modelcontainer',
            name='fully_trained',
            field=models.ForeignKey(
                blank=True, to='DRP.StatsModel', null=True),
        ),
        migrations.AlterField(
            model_name='modelcontainer',
            name='modelVisitorLibrary',
            field=models.CharField(max_length=200),
        ),
        migrations.AlterField(
            model_name='modelcontainer',
            name='modelVisitorTool',
            field=models.CharField(max_length=200),
        ),
        migrations.AlterField(
            model_name='modelcontainer',
            name='splitter',
            field=models.CharField(default=b'', max_length=200, blank=True),
        ),
    ]
