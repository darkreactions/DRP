'''A module containing views pertinent to the manipulation of reaction objects'''

from django.views.generic import CreateView, ListView, UpdateView
from DRP.models import PerformedReaction
from django.utils.decorators import method_decorator
from decorators import userHasLabGroup, hasSignedLicense, labGroupSelected
from django.contrib.auth.decorators import login_required

class ListPerformedReactions(ListView):

  template_name='reactions_list.html'
  context_object_name='reactions'
  model=PerformedReaction 
  
  @method_decorator(login_required)
  @method_decorator(hasSignedLicense)
  @method_decorator(userHasLabgroup)
  @labGroupSelected #sets self.labGroup
  def dispatch(self, request, *args, **kwargs)
    self.queryset = self.labGroup.performedreaction_set.all()
    return super(ListCompound, self).dispatch(request, *args, **kwargs)

  def get_context_data(self, **kwargs):
    context = super(ListCompound, self).get_context_data(**kwargs)
    context['lab_form'] = self.labForm
#    context['filter_formset'] = self.filterFormSet
    return context
