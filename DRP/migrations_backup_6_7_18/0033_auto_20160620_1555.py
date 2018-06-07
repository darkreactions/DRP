# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.db.models import Count


def rearrange_compound_guides(apps, schema_editor):
    """
    Rearrange the data for the compound guides so that abbreviations and compounds have a layer of indirection.

    This permits us a better storage of compounds and their relationships with groups in an unambigious way.
    """
    CompoundQuantity = apps.get_model('DRP', "CompoundQuantity")
    Compound = apps.get_model("DRP", "Compound")
    CompoundGuideEntry = apps.get_model("DRP", "CompoundGuideEntry")
    dups = Compound.objects.exclude(CSID=None).values('CSID').annotate(
        csid_count=Count('id')).order_by().filter(csid_count__gt=1)
    for compound in Compound.objects.all():
        if compound.CSID is not None:
            actualCompound = Compound.objects.filter(CSID=compound.CSID)[
                0]  # Ignore duplicate compounds
        else:
            # We don't want to accidentally strip out all of the custom
            # compounds.
            acutalCompound = compound
        if not CompoundGuideEntry.objects.filter(compound=actualCompound).exists():
            CompoundGuideEntry.objects.create(
                compound=actualCompound, labGroup=compound.labGroup, abbrev=compound.abbrev)
    for csidValue in dups:  # This code sets all of the compound quantities to use the same set of compound
        csid = csidValue['CSID']
        compounds = Compound.objects.filter(CSID=csid)
        for compoundQuantity in CompoundQuantity.objects.filter(compound__CSID=csid):
            compoundQuantity.compound = compounds[0]
            compoundQuantity.save()
        for compound in compounds[1:]:  # Delete the known duplicates.
            compound.delete()


class Migration(migrations.Migration):

    dependencies = [
        ('DRP', '0032_auto_20160620_1554'),
    ]

    operations = [
        migrations.RunPython(rearrange_compound_guides)
    ]
