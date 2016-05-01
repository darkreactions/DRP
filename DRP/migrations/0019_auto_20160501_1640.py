# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0018_auto_20160501_1633'),
    ]

    operations = [
        migrations.AlterField(
            model_name='performedreaction',
            name='performedBy',
            field=models.ForeignKey(related_name='performedReactions', default=None, blank=True, to=settings.AUTH_USER_MODEL, null=True),
        ),
        migrations.AlterField(
            model_name='performedreaction',
            name='performedDateTime',
            field=models.DateTimeField(default=None, help_text=b'Date in format YYYY-MM-DD', null=True, verbose_name=b'Date Reaction Performed', blank=True),
        ),
    ]
