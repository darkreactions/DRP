'''A module containing views pertinent to compound objects'''

from django.contrib.auth.models import User
from django.views.generic import CreateView, ListView
from DRP.models import Compound
from DRP.forms import CompoundForm
from django.utils.decorators import method_decorator
from decorators import userHasLabGroup


class CreateCompound(CreateView):
  '''A view managing the creation of compound objects'''

  model=Compound
  form_class = CompoundForm
  template_name='compound_form.html'
 
  def get_form_kwargs(self):
    '''Overridden to add the request.user value into the kwargs'''
    kwargs = super(CreateCompound, self).get_form_kwargs() 
    kwargs['user']=self.request.user
    return kwargs

  @method_decorator(userHasLabGroup)
  def dispatch(self, request, *args, **kwargs):
    '''Overridden with a decorator to ensure that a user is at least logged in'''
    return super(CreateCompound, self).dispatch(request, *args, **kwargs)

  def get_context_data(self, **kwargs):
    context = super(CreateCompound, self).get_context_data(**kwargs)
    context['page_heading'] = 'Add a New Compound'
    return context
    
class ListCompound(ListView):
  '''A view managing the viewing of the compound guide'''

  template_name='compound_list.html'
  context_object_name='compounds'

  @method_decorator(userHasLabGroup)
  def dispatch(self, request, *args, **kwargs):
    '''Overriden with a decorator to ensure that user is logged in and has at least one labGroup
    Relates the queryset of this view to the logged in user.
    '''

    if request.user.labgroup_set.all().count() > 1:
      if 'labgroup_id' in request.session and request.user.labgroup_set.filter(pk=request.session['labgroup_id']).count() > 1:
        self.queryset = request.user.labgroup_set.get(pk=request_session['labgroup_id']).compound_set.all()
      else:
        self.queryset = []
    else:
      #user only has one labgroup, so don't bother asking which group's compoundlist they want to look at.
      self.queryset = request.user.labgroup_set.all()[0].compound_set.all()
    return super(ListCompound, self).dispatch(request, *args, **kwargs)
