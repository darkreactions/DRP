# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0022_auto_20160510_0941'),
    ]

    operations = [
        migrations.AlterField(
            model_name='descriptor',
            name='name',
            field=models.CharField(max_length=300, verbose_name=b'Full name'),
        ),
    ]
