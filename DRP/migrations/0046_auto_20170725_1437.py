# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.performedReaction
import DRP.models.fileStorage


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0045_auto_20160928_0849'),
    ]

    operations = [
        migrations.AlterField(
            model_name='performedreaction',
            name='labBookPage',
            field=models.ImageField(null=True, storage=DRP.models.fileStorage.OverwriteStorage(location='/home/h205c/DRP/sec_media/lab_notes', base_url='/database/lab_notes/'), verbose_name='Lab Book Image', upload_to=DRP.models.performedReaction.UploadLabNotesTo()),
        ),
    ]
