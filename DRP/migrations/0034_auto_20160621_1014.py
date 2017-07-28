# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0033_auto_20160620_1555'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='compound',
            unique_together=set([('CSID', 'labGroup')]),
        ),
    ]
