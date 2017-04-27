"""A module containing urls for the database (reactions, compound guide) components of DRP."""

from django.conf.urls import url
from django.views.defaults import permission_denied
import DRP.views
from DRP.models import NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue
from DRP.models import OrdRxnDescriptorValue
from DRP.forms import NumRxnDescValFormFactory, OrdRxnDescValFormFactory, BoolRxnDescValFormFactory, CatRxnDescValFormFactory

urls = [
    url('^(?P<filetype>.csv|.html|.arff)?$',
        DRP.views.reaction.ListPerformedReactions.as_view(), name='reactionlist_typed'),
    url('^/$', DRP.views.reaction.ListPerformedReactions.as_view(), name='reactionlist'),
    url('^/add.html', DRP.views.reaction.createReaction, name='newReaction'),
    url('^/entry_(?P<rxn_id>\d+)/compoundquantities.html',
        DRP.views.reaction.addCompoundDetails, name="addCompoundDetails"),
    url('^/entry_(?P<rxn_id>\d+)/num_desc_vals.html', DRP.views.reaction.createGenDescVal,
        {'descValClass': NumRxnDescriptorValue,
         'descValFormClass': NumRxnDescValFormFactory,
         'infoHeader': 'Numerical Descriptor Values',
         'createNext': 'createOrdDescVals'
         },
        name="createNumDescVals"),
    url('^/entry_(?P<rxn_id>\d+)/ord_desc_vals.html', DRP.views.reaction.createGenDescVal,
        {'descValClass': OrdRxnDescriptorValue,
         'descValFormClass': OrdRxnDescValFormFactory,
         'infoHeader': 'Ordinal Descriptor Values',
         'createNext': 'createBoolDescVals'
         },
        name="createOrdDescVals"),
    url('^/entry_(?P<rxn_id>\d+)/cat_desc_vals.html', DRP.views.reaction.createGenDescVal,
        {'descValClass': CatRxnDescriptorValue,
         'descValFormClass': CatRxnDescValFormFactory,
         'infoHeader': 'Categorical Descriptor Values',
         'createNext': None
         },
        name="createCatDescVals"),
    url('^/entry_(?P<rxn_id>\d+)/bool_desc_vals.html', DRP.views.reaction.createGenDescVal,
        {'descValClass': BoolRxnDescriptorValue,
         'descValFormClass': BoolRxnDescValFormFactory,
         'infoHeader': 'Boolean Descriptor Values',
         'createNext': 'createCatDescVals'
         },
        name="createBoolDescVals"),
    url('^/entry_(?P<rxn_id>\d+)/display', DRP.views.reaction.displayReaction, name='reactionInfo'),
    url('^/entry_(?P<rxn_id>\d+)/',
        DRP.views.reaction.editReaction, name='editReaction'),
    url('^/delete$', DRP.views.reaction.deleteReaction, name='deleteReaction'),
    url('^/invalidate$', DRP.views.reaction.invalidateReaction,
        name='invalidateReaction'),
    url('^/lab_notes/$', permission_denied),
    url('^/lab_notes/(?P<labgroup_id>\d+)_(?P<reference>[^//]*).jpg',
        DRP.views.reaction.labBookImage, name="labBookImage"),
    url('^/import/apiv1/(?P<component>[^//]*).xml', DRP.views.api1),
    url('^/select_viewing_group.html', DRP.views.selectGroup, name='selectGroup'),
    url('^/compoundguide(?P<filetype>.csv|.html|.arff|/)$',
        DRP.views.compound.ListCompound.as_view(), name='compoundguide'),
    url('^/compoundguide/add.html$',
        DRP.views.compound.CreateCompound.as_view(), name='newCompound'),
    url('^/compoundguide/delete$',
        DRP.views.compound.deleteCompound, name='deleteCompound'),
    url('^/compoundguide/edit_(?P<pk>\d+).html',
        DRP.views.compound.EditCompound.as_view(), name='editCompound'),
    url('^/compoundguide/display_(?P<pk>\d+).html',
        DRP.views.compound.displayCompound, name='displayCompound')
]
