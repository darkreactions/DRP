# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0044_auto_20160928_0835'),
    ]

    operations = [
        migrations.RunSQL(
            'CREATE TRIGGER soil_compound BEFORE UPDATE ON DRP_compound FOR EACH ROW SET NEW.dirty = 1, NEW.recalculate = OLD.calculating;',
            'DROP TRIGGER soil_compound;'
        ),
        migrations.RunSQL(
            'CREATE TRIGGER soil_reaction BEFORE UPDATE ON DRP_reaction FOR EACH ROW SET NEW.dirty = 1, NEW.recalculate = OLD.calculating;',
            'DROP TRIGGER soil_reaction;'
        )
    ]
