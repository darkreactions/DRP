# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations

def fix_sql_mess(apps, schema_editor):
    Compound = apps.get_model("DRP", "Compound") 
    constraints = schema_editor._constraint_names(Compound, ["labGroup_id", "CSID",], unique=True)
    if len(constraints) > 0:
        schema_editor.execute(schema_editor._delete_constraint_sql(
                schema_editor.sql_delete_index,
                Compound,
                constraints[0]
            ))


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0034_auto_20160621_1014'),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            [migrations.RunPython(fix_sql_mess)],
            [migrations.AlterUniqueTogether(name='compound', unique_together=set([]))]
        )
    ]
