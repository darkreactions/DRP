# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.performedReaction
import django.core.files.storage
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0037_performedreaction_labbookpage'),
    ]

    operations = [
        migrations.AlterField(
            model_name='performedreaction',
            name='duplicateOf',
            field=models.ForeignKey(blank=True, verbose_name='Duplicate Of', default=None,
                                    related_name='duplicatedBy', to='DRP.PerformedReaction', null=True),
        ),
        migrations.AlterField(
            model_name='performedreaction',
            name='labBookPage',
            field=models.ImageField(upload_to=DRP.models.performedReaction.UploadLabNotesTo(), verbose_name='Lab Book Image', storage=django.core.files.storage.FileSystemStorage(
                base_url='/database/lab_notes/', location='/home/padler1/programming/drp/DRP/sec_media/lab_notes'), null=True),
        ),
        migrations.AlterField(
            model_name='performedreaction',
            name='performedBy',
            field=models.ForeignKey(blank=True, verbose_name='Performed By', default=None,
                                    related_name='performedReactions', to=settings.AUTH_USER_MODEL, null=True),
        ),
    ]
