'''A module containing views pertinent to compound objects'''

from django.contrib.auth.models import User
from django.views.generic import CreateView
from DRP.models import Compound
from DRP.forms import CompoundForm
from django.contrib.auth.decorators import login_required
from django.utils.decorators import method_decorator

class CreateCompound(CreateView):

  model=Compound
  form_class = CompoundForm
  template_name='compound_form.html'
 
  def get_form_kwargs(self):
    '''Overridden to add the request.user value into the kwargs'''
    kwargs = super(CreateCompound, self).get_form_kwargs() 
    kwargs['user']=self.request.user
    return kwargs

  @method_decorator(login_required)
  def dispatch(self, *args, **kwargs):
    '''Overridden with a decorator to ensure that a user is at least logged in'''
    return super(CreateCompound, self).dispatch(*args, **kwargs)
    
