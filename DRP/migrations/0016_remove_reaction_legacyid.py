# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0015_reaction_legacyid'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='reaction',
            name='legacyID',
        ),
    ]
