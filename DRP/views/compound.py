'''A module containing views pertinent to compound objects'''

from django.contrib.auth.models import User
from django.views.generic import CreateView, ListView, UpdateView
from DRP.models import Compound
from DRP.forms import CompoundForm, LabGroupSelectionForm, CompoundEditForm, CompoundDeleteForm
from django.utils.decorators import method_decorator
from decorators import userHasLabGroup, hasSignedLicense
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse_lazy as reverse
from django.shortcuts import redirect
from django.utils.http import urlencode
from django.http import HttpResponse, Http404, HttpResponseForbidden
from django.template.loader import get_template
from django.views.decorators.http import require_POST
from django.template import RequestContext

class CreateCompound(CreateView):
  '''A view managing the creation of compound objects'''

  model=Compound
  form_class = CompoundForm
  template_name='compound_form.html'
  success_url=reverse('compoundguide')
 
  def get_form_kwargs(self):
    '''Overridden to add the request.user value into the kwargs'''
    kwargs = super(CreateCompound, self).get_form_kwargs() 
    kwargs['user']=self.request.user
    return kwargs

  @method_decorator(login_required)
  @method_decorator(hasSignedLicense)
  @method_decorator(userHasLabGroup)
  def dispatch(self, request, *args, **kwargs):
    '''Overridden with a decorator to ensure that a user is at least logged in'''
    return super(CreateCompound, self).dispatch(request, *args, **kwargs)

  def get_context_data(self, **kwargs):
    context = super(CreateCompound, self).get_context_data(**kwargs)
    context['page_heading'] = 'Add a New Compound'
    return context
    
class EditCompound(UpdateView):
  '''A view managing the editing of compound objects'''

  form_class=CompoundEditForm
  template_name='compound_edit.html'
  success_url=reverse('compoundguide')
  model = Compound

  @method_decorator(login_required)
  @method_decorator(hasSignedLicense)
  @method_decorator(userHasLabGroup)
  def dispatch(self, request, *args, **kwargs):
    '''Checks user has sufficient credentials and has row-level permissions for this compound'''
    try:
      compound = Compound.objects.get(pk=self.get_object().pk, labGroup__in=request.user.labgroup_set.all())
    except Compound.DoesNotExist:
      raise Http404("A compound matching your query could not be found.")
    if compound.custom:
      template = get_template('compound_403.html')
      return HttpResponseForbidden(template.render(RequestContext(request)))
    else:
      return super(EditCompound, self).dispatch(request, *args, **kwargs)

@require_POST
@login_required
@hasSignedLicense
@userHasLabGroup
def deleteCompound(request, *args, **kwargs):
  '''A view managing the deletion of compound objects'''
  form = CompoundDeleteForm(data=request.POST, user=request.user) 
  if form.is_valid():
    form.save()
  return redirect('compoundguide')

class ListCompound(ListView):
  '''A view managing the viewing of the compound guide'''

  template_name='compound_list.html'
  context_object_name='compounds'
  model=Compound

  @method_decorator(login_required)
  @method_decorator(hasSignedLicense)
  @method_decorator(userHasLabGroup)
  def dispatch(self, request, *args, **kwargs):
    '''Overriden with a decorator to ensure that user is logged in and has at least one labGroup
    Relates the queryset of this view to the logged in user.
    '''
    
    self.lab_form = LabGroupSelectionForm(request.user)
    if request.user.labgroup_set.all().count() > 1:
      if 'labgroup_id' in request.session and request.user.labgroup_set.filter(pk=request.session['labgroup_id']).exists():
        self.queryset = request.user.labgroup_set.get(pk=request.session['labgroup_id']).compound_set.all()
        self.lab_form.fields['labGroup'].initial = request.session['labgroup_id']
        return super(ListCompound, self).dispatch(request, *args, **kwargs)
      elif 'labgroup_id' not in request.session:
        return redirect(reverse('selectGroup') + '?{0}'.format(urlencode({'next':request.path_info})))
      else:
        raise RuntimeError("This shouldn't happen")
    else:
      #user only has one labgroup, so don't bother asking which group's compoundlist they want to look at.
      self.queryset = request.user.labgroup_set.all()[0].compound_set.all()
      return super(ListCompound, self).dispatch(request, *args, **kwargs)

  def get_context_data(self, **kwargs):
    context = super(ListCompound, self).get_context_data(**kwargs)
    context['lab_form'] = self.lab_form
    return context

