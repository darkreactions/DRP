# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.validators


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0025_auto_20160524_1613'),
    ]

    operations = [
        migrations.AlterField(
            model_name='performedreaction',
            name='performedDateTime',
            field=models.DateTimeField(default=None, validators=[DRP.models.validators.notInTheFuture], blank=True,
                                       help_text='Timezone assumed EST, Date in format YYYY-MM-DD', null=True, verbose_name='Date Reaction Performed'),
        ),
    ]
