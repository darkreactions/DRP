# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0020_auto_20160501_1748'),
    ]

    operations = [
        migrations.AlterField(
            model_name='performedreaction',
            name='convertedLegacyRef',
            field=models.CharField(blank=True, max_length=40, null=True, validators=[django.core.validators.RegexValidator(
                b'[a-z0-9._]*[a-z][a-z0-9._]*', b'Please include only values which are limited to alphanumeric characters, underscores, periods, and must include at least one alphabetic character.')]),
        ),
        migrations.AlterField(
            model_name='performedreaction',
            name='reference',
            field=models.CharField(max_length=40, validators=[django.core.validators.RegexValidator(
                b'[a-z0-9._]*[a-z][a-z0-9._]*', b'Please include only values which are limited to alphanumeric characters, underscores, periods, and must include at least one alphabetic character.')]),
        ),
    ]
