'''A module containing urls for the database (reactions, compound guide) components of DRP'''

from django.conf.urls import patterns, include, url
from DRP.models import Compound
import DRP.views

urls = patterns('',
  url('^select_viewing_group.html', DRP.views.selectGroup, name='selectGroup'),
  url('^compoundguide(?P<filetype>.csv|.html|/)$', DRP.views.compound.ListCompound.as_view(), name='compoundguide'),
  url('^compoundguide/search(?P<filetype>.html|.csv)$', DRP.views.compound.ListCompound.as_view(), name='compoundSearch'),
  url('^compoundguide/advanced_search(?P<filetype>.html|.csv)$', DRP.views.compound.AdvancedCompoundSearchView.as_view(), name='advCompoundSearch'),
  url('^compoundguide/add.html$', DRP.views.compound.CreateCompound.as_view(), name='newCompound'),
  url('^compoundguide/delete$', DRP.views.compound.deleteCompound, name='deleteCompound'),
  url('^compoundguide/edit_(?P<pk>\d+).html', DRP.views.compound.EditCompound.as_view(), name='editCompound'),
  url('^compoundguide/upload.html', DRP.views.compound.uploadCompound, name='uploadcompoundcsv')
)
