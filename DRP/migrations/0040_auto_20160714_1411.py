# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.performedReaction
import DRP.models.fileStorage


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0039_auto_20160630_1321'),
    ]

    operations = [
        migrations.AlterField(
            model_name='performedreaction',
            name='labBookPage',
            field=models.ImageField(verbose_name='Lab Book Image', null=True, upload_to=DRP.models.performedReaction.UploadLabNotesTo(), storage=DRP.models.fileStorage.OverwriteStorage(location='/home/mdbyars/DRP/sec_media/lab_notes', base_url='/database/lab_notes/')),
        ),
    ]
