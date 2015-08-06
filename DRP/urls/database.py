'''A module containing urls for the database (reactions, compound guide) components of DRP'''

from django.conf.urls import patterns, include, url
from django.views.generic.edit import CreateView
from DRP.models import Compound
from DRP.forms import CompoundForm
from django.contrib.auth.decorators import login_required
import DRP.views

urls = patterns('',
#  url('^compoundguide/$', DRP.views.compoundguide.view),
  url('^compoundguide/add.html$',
      login_required(CreateView.as_view(model=Compound,
      template_name='compound_form.html',
      form_class=CompoundForm)),
      name='newCompound')
#  url('^compoundguide/delete$', DRP.views.compoundguide.delete),
#  url('^compoundguide/(?P<compound_id>\d+).html'), DRP.views.compoundguide.edit)
)
