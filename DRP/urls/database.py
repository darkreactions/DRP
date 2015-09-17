'''A module containing urls for the database (reactions, compound guide) components of DRP'''

from django.conf.urls import patterns, include, url
import DRP.views

urls = patterns('',
  url('^$', DRP.views.reaction.ListPerformedReactions.as_view(), name='reactionlist'),
  url('^add.html', DRP.views.reaction.reactionForm, name='newReaction'),
  url('^edit_(?P<pk>\d+).html', DRP.views.reaction.reactionForm, name='editReaction'),
  url('^select_viewing_group.html', DRP.views.selectGroup, name='selectGroup'),
  url('^compoundguide(?P<filetype>.csv|.html|.arff|/)$', DRP.views.compound.ListCompound.as_view(), name='compoundguide'),
  url('^compoundguide/search(?P<filetype>.html|.csv|.arff)$', DRP.views.compound.ListCompound.as_view(), name='compoundSearch'),
  url('^compoundguide/advanced_search(?P<filetype>.html|.csv|.arff)$', DRP.views.compound.AdvancedCompoundSearchView.as_view(), name='advCompoundSearch'),
  url('^compoundguide/add.html$', DRP.views.compound.CreateCompound.as_view(), name='newCompound'),
  url('^compoundguide/delete$', DRP.views.compound.deleteCompound, name='deleteCompound'),
  url('^compoundguide/edit_(?P<pk>\d+).html', DRP.views.compound.EditCompound.as_view(), name='editCompound'),
  url('^compoundguide/upload.html', DRP.views.compound.uploadCompound, name='uploadcompoundcsv')
)
