# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0024_auto_20160512_1658'),
    ]

    operations = [
        migrations.AlterField(
            model_name='performedreaction',
            name='convertedLegacyRef',
            field=models.CharField(blank=True, max_length=40, null=True, validators=[django.core.validators.RegexValidator(
                '^[a-z0-9._]*[a-z][a-z0-9._]*$', 'Please include only values which are limited to alphanumeric characters, underscores, periods, and must include at least one alphabetic character.')]),
        ),
        migrations.AlterField(
            model_name='performedreaction',
            name='performedDateTime',
            field=models.DateTimeField(default=None, help_text='Timezone assumed EST, Date in format YYYY-MM-DD',
                                       null=True, verbose_name='Date Reaction Performed', blank=True),
        ),
        migrations.AlterField(
            model_name='performedreaction',
            name='reference',
            field=models.CharField(max_length=40, validators=[django.core.validators.RegexValidator(
                '^[a-z0-9\\._]*[a-z][a-z0-9\\._]*$', 'Please include only values which are limited to alphanumeric characters, underscores, periods, and must include at least one alphabetic character.')]),
        ),
    ]
