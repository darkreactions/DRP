# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0028_auto_20160526_1446'),
    ]

    operations = [
        migrations.AlterField(
            model_name='boolmoldescriptorvalue',
            name='value',
            field=models.NullBooleanField(
                verbose_name=b'Value for descriptor'),
        ),
        migrations.AlterField(
            model_name='boolrxndescriptorvalue',
            name='value',
            field=models.NullBooleanField(
                verbose_name=b'Value for descriptor'),
        ),
        migrations.AlterField(
            model_name='catmoldescriptorvalue',
            name='value',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT,
                                    blank=True, to='DRP.CategoricalDescriptorPermittedValue', null=True),
        ),
        migrations.AlterField(
            model_name='catrxndescriptorvalue',
            name='value',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT,
                                    blank=True, to='DRP.CategoricalDescriptorPermittedValue', null=True),
        ),
        migrations.AlterField(
            model_name='nummoldescriptorvalue',
            name='value',
            field=models.FloatField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='numrxndescriptorvalue',
            name='value',
            field=models.FloatField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='ordmoldescriptorvalue',
            name='value',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='ordrxndescriptorvalue',
            name='value',
            field=models.IntegerField(null=True, blank=True),
        ),
    ]
