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
        migrations.AddField(
            model_name='labgroup',
            name='defaultDescriptors',
            field=models.ManyToManyField(related_name='isDefaultForLabGroups', to='DRP.Descriptor'),
        ),
        migrations.AlterField(
            model_name='performedreaction',
            name='labBookPage',
            field=models.ImageField(null=True, storage=DRP.models.fileStorage.OverwriteStorage(base_url='/database/lab_notes/', location='/home/padler1/programming/drp/DRP/sec_media/lab_notes'), upload_to=DRP.models.performedReaction.UploadLabNotesTo(), verbose_name='Lab Book Image'),
        ),
    ]
