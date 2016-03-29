# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0006_auto_20160302_2159'),
    ]

    operations = [
        #migrations.DeleteModel(
        #    name='CompoundQuantityManager',
        #),
        migrations.RemoveField(
            model_name='modelcontainer',
            name='featureLibrary',
        ),
        migrations.RemoveField(
            model_name='modelcontainer',
            name='featureTool',
        ),
        migrations.AlterField(
            model_name='descriptor',
            name='calculatorSoftware',
            field=models.CharField(max_length=100, validators=[django.core.validators.RegexValidator(b'[A-Za-z0-9][A-Za-z0-9_]+', b'Please include only values which are limited toalphanumeric characters and underscores, and must startwith an alphabetic character.')]),
        ),
        migrations.AlterField(
            model_name='descriptor',
            name='calculatorSoftwareVersion',
            field=models.CharField(max_length=20, validators=[django.core.validators.RegexValidator(b'[A-Za-z0-9][A-Za-z0-9_]+', b'Please include only values which are limited toalphanumeric characters and underscores, and must startwith an alphabetic character.')]),
        ),
    ]
