# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0010_auto_20160330_2126'),
    ]

    operations = [
        migrations.AddField(
            model_name='statsmodel',
            name='inputFile',
            field=models.FileField(
                max_length=255, upload_to='model_inputs', blank=True),
        ),
    ]
