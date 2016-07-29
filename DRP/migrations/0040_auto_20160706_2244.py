# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.fileStorage
import DRP.models.performedReaction


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0039_auto_20160630_1321'),
    ]

    operations = [
        migrations.AlterField(
            model_name='performedreaction',
            name='labBookPage',
            field=models.ImageField(storage=DRP.models.fileStorage.OverwriteStorage(base_url='/database/lab_notes/', location='/home/rich/DRP/sec_media/lab_notes'), verbose_name='Lab Book Image', upload_to=DRP.models.performedReaction.UploadLabNotesTo(), null=True),
        ),
    ]
