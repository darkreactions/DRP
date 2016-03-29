# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0009_auto_20160329_1826'),
    ]

    operations = [
        migrations.AlterField(
            model_name='descriptor',
            name='calculatorSoftware',
            field=models.CharField(max_length=100, validators=[django.core.validators.RegexValidator(b'[A-Za-z0-9][A-Za-z0-9_]*', b'Please include only values which are limited to alphanumeric characters and underscores, and must start with an alphabetic character.')]),
        ),
        migrations.AlterField(
            model_name='descriptor',
            name='calculatorSoftwareVersion',
            field=models.CharField(max_length=20, validators=[django.core.validators.RegexValidator(b'[A-Za-z0-9][A-Za-z0-9_.]*', b'Please include only values which are limited to alphanumeric characters, periods and underscores, and must start with an alphabetic character.')]),
        ),
        migrations.AlterField(
            model_name='descriptor',
            name='heading',
            field=models.CharField(max_length=200, validators=[django.core.validators.RegexValidator(b'[A-Za-z0-9][A-Za-z0-9_]*', b'Please include only values which are limited to alphanumeric characters and underscores, and must start with an alphabetic character.')]),
        ),
    ]
