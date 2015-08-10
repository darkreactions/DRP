'''A module containing urls for the database (reactions, compound guide) components of DRP'''

from django.conf.urls import patterns, include, url
from DRP.models import Compound
import DRP.views

urls = patterns('',
  url('^select_viewing_group.html', DRP.views.database.selectGroup, name='selectGroup'),
  url('^compoundguide/$', DRP.views.compound.ListCompound.as_view(), name='compoundguide'),
  url('^compoundguide/add.html$', DRP.views.compound.CreateCompound.as_view(), name='newCompound')
#  url('^compoundguide/delete$', DRP.views.compoundguide.delete),
#  url('^compoundguide/(?P<compound_id>\d+).html'), DRP.views.compoundguide.edit)
)
