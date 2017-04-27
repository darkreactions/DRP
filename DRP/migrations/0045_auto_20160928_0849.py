# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0044_auto_20160928_0835'),
    ]

    operations = [
        migrations.RunSQL(
            'CREATE TRIGGER soil_compound BEFORE UPDATE ON DRP_compound FOR EACH ROW SET NEW.dirty = IF(OLD.calculating=1 and NEW.calculating=0, 0, 1), NEW.recalculate = OLD.calculating AND NEW.calculating AND (OLD.recalculate AND NEW.recalculate);',
            'DROP TRIGGER soil_compound;'
        ),
        migrations.RunSQL(
            'CREATE TRIGGER soil_compound_reaction AFTER UPDATE ON DRP_compound FOR EACH ROW UPDATE DRP_reaction JOIN DRP_compoundquantity ON DRP_reaction.id=DRP_compoundquantity.reaction_id SET DRP_reaction.dirty = 1 WHERE DRP_compoundquantity.compound_id=NEW.id;',
            'DROP TRIGGER soil_compound_reaction;'
        ),
        migrations.RunSQL(
            'CREATE TRIGGER soil_reaction BEFORE UPDATE ON DRP_reaction FOR EACH ROW SET NEW.dirty = IF(OLD.calculating=1 and NEW.calculating=0, 0, 1), NEW.recalculate = OLD.calculating AND NEW.calculating AND (OLD.recalculate AND NEW.recalculate);',
            'DROP TRIGGER soil_reaction;'
        )
    ]
