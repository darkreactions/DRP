'''A module containing views pertinent to the manipulation of reaction objects'''

from django.views.generic import CreateView, ListView, UpdateView
from DRP.models import PerformedReaction, CompoundQuantity
from DRP.forms import PerformedRxnForm
from DRP.forms import NumRxnDescValForm, OrdRxnDescValForm, BoolRxnDescValForm, CatRxnDescValForm  
from django.utils.decorators import method_decorator
from decorators import userHasLabGroup, hasSignedLicense, labGroupSelected
from django.contrib.auth.decorators import login_required
from django.forms.models import modelformset_factory
from DRP.forms import ModelFormSet, FormSet
from django.forms.formsets import TOTAL_FORM_COUNT
from django.shortcuts import render, redirect

class ListPerformedReactions(ListView):

  template_name='reactions_list.html'
  context_object_name='reactions'
  model=PerformedReaction 
  
  @method_decorator(login_required)
  @method_decorator(hasSignedLicense)
  @method_decorator(userHasLabGroup)
  @labGroupSelected #sets self.labGroup
  def dispatch(self, request, *args, **kwargs):
    self.queryset = PerformedReaction.objects.filter(reaction_ptr__in=self.labGroup.reaction_set.all()) | PerformedReaction.objects.filter(public=True)
    return super(ListPerformedReactions, self).dispatch(request, *args, **kwargs)

  def get_context_data(self, **kwargs):
    context = super(ListPerformedReactions, self).get_context_data(**kwargs)
    context['lab_form'] = self.labForm
#    context['filter_formset'] = self.filterFormSet
    return context

@login_required
@hasSignedLicense
@userHasLabGroup
def createReaction(request):
  descFields = ('descriptor', 'value')
  if request.method=='POST':
    reactantsFormSetInst = ModelFormSet(CompoundQuantity, fields=('compound', 'role', 'amount'), data=request.POST, canAdd=True, canDelete=True)
    reactionForm = PerformedRxnForm(request.user, data=request.POST) 

    descriptorFormSets = (
      FormSet(NumRxnDescValForm, data=request.POST, prefix='num', canDelete=True),
      FormSet(OrdRxnDescValForm, data=request.POST, prefix='ord', canDelete=True),
      FormSet(CatRxnDescValForm, data=request.POST, prefix='cat', canDelete=True),
      FormSet(BoolRxnDescValForm, data=request.POST, prefix='bool', canDelete=True)
    )

    if 'save' in request.POST:
      if reactionForm.is_valid() and reactantsFormSetInst.is_valid() and all(d.is_valid() for d in descriptorFormSets):
        rxn = reactionForm.save()
        for reactants in reactantsFormSetInst.save(commit=False):
          reactant.reaction=rxn.reaction
          reactant.save()
        for formSet in descriptorFormSets:
          for descriptorValue in formset.save(commit=False):
            descriptorValue.reaction=rxn.reaction
        return redirect('reactionlist')
  else:
    reactionForm = PerformedRxnForm(request.user)
    reactantsFormSetInst = ModelFormSet(CompoundQuantity, fields=('compound', 'role', 'amount'), canAdd=True)
    descriptorFormSets = (
      FormSet(NumRxnDescValForm, prefix='num'),
      FormSet(OrdRxnDescValForm,prefix='ord'),
      FormSet(BoolRxnDescValForm, prefix='bool'),
      FormSet(CatRxnDescValForm, prefix='cat')
    )
  return render(request, 'reaction_form.html', {'reaction_form':reactionForm, 'reactants_formset':reactantsFormSetInst, 'descriptor_formsets':descriptorFormSets}) 
