# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0031_auto_20160620_1357'),
    ]

    operations = [
        migrations.CreateModel(
            name='CompoundGuideEntry',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('abbrev', models.CharField(
                    max_length=100, verbose_name='Abbreviation')),
                ('compound', models.ForeignKey(to='DRP.Compound')),
                ('labGroup', models.ForeignKey(to='DRP.LabGroup')),
            ],
        ),
        migrations.AddField(
            model_name='compound',
            name='labGroups',
            field=models.ManyToManyField(
                to='DRP.LabGroup', verbose_name='Lab Groups', through='DRP.CompoundGuideEntry'),
        ),
        migrations.AlterUniqueTogether(
            name='compoundguideentry',
            unique_together=set(
                [('compound', 'labGroup'), ('abbrev', 'labGroup')]),
        ),
    ]
