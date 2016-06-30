# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.core.files.storage
import DRP.models.performedReaction


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0036_auto_20160621_1220'),
    ]

    operations = [
        migrations.AddField(
            model_name='performedreaction',
            name='labBookPage',
            field=models.ImageField(null=True, upload_to=DRP.models.performedReaction.UploadLabNotesTo(), storage=django.core.files.storage.FileSystemStorage(base_url='/database/lab_notes/', location='/home/padler1/programming/drp/DRP/sec_media/lab_notes')),
        ),
    ]
