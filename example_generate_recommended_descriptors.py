import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'
import django
# import time

django.setup()

from django.core.management.base import BaseCommand
from DRP.models import Reaction, Compound, NumRxnDescriptorValue
from django import db
from django.conf import settings
import logging
import importlib
from django.db import transaction

import DRP

rxnDescriptorPlugins = [importlib.import_module(plugin) for
                        plugin in settings.RXN_DESCRIPTOR_PLUGINS]



underpopulated_pks = [i.pk for i in Reaction.objects.all() if NumRxnDescriptorValue.objects.filter(reaction__pk=i.pk).count() >= 500 ]

#reactions = DRP.models.PerformedReaction.objects.filter(pk__in=underpopulated_pks)





for i in range(30):
  reactions = Reaction.objects.filter(pk__in=(underpopulated_pks+underpopulated_pks)[i*len(underpopulated_pks):i+1*len(underpopulated_pks)])
  for plugin in rxnDescriptorPlugins[1:]:
      print(plugin)
      plugin.calculate_many(reactions, whitelist=csvheaders, bulk_delete=True)




# print(DRP.models.NumRxnDescriptorValue.objects.filter(reaction__pk=36449).count())
# print('finished ' + str(plugin) + ' at :' + str(time.time()))
#
#
# not_c = []
# for rec_desc_header in rec_descs:
#     if DRP.models.NumRxnDescriptorValue.objects.filter(reaction__pk=36449, descriptor__heading=rec_desc_header).count() == 1:
#         pass
#     elif DRP.models.BoolRxnDescriptorValue.objects.filter(reaction__pk=36449, descriptor__heading=rec_desc_header).count() == 1:
#         pass
#     elif DRP.models.OrdRxnDescriptorValue.objects.filter(reaction__pk=36449, descriptor__heading=rec_desc_header).count() == 1:
#         pass
#     else:
#         not_c.append(rec_desc_header)
#
# print(not_c)
#
