'''A module containing views pertinent to the manipulation of reaction objects'''

from django.views.generic import CreateView, ListView, UpdateView
from DRP.models import PerformedReaction, CompoundQuantity
from DRP.forms import PerformedRxnForm
from DRP.forms import NumRxnDescValFormSet, OrdRxnDescValFormSet, BoolRxnDescValFormSet, CatRxnDescValFormSet  
from django.utils.decorators import method_decorator
from decorators import userHasLabGroup, hasSignedLicense, labGroupSelected
from django.contrib.auth.decorators import login_required
from django.forms import modelformset_factory
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
    reactionForm = PerformedRxnForm(request.user, data=request.POST) 

    reactantsFormSet = modelformset_factory(CompoundQuantity, fields=('compound', 'role', 'amount'))(data=request.POST)
    if 'add_reactant' in request.POST:
      reactantsFormSet.extra = reactantsFormSet.total_form_count()+1
    if len(reactantsFormset.forms) > 1:
      reactantsFormSet.can_delete=True

    descriptorFormsets = (
      NumRxndescValFormSet(data=request.POST), OrdRxndescValFormSet(data=request.POST), BoolRxnDescValFormSet(data=request.POST), CatRxnDescValFormSet(data=request.POST)
    )
    for formSet in descriptorFormSets:
      if 'add_' + formSet.prefix + '_descriptor' in request.POST:
        formSet.extra = reactantsFormSet.total_form_count()+1
    if reactionForm.is_valid() and reactantsFormSet.is_valid() and all(d.is_valid() for d in descriptorFormsets):
      rxn = reactionForm.save()
      for reactants in reactantsFormSet.save(commit=False):
        reactant.reaction=rxn.reaction
        reactant.save()
      for formSet in descriptorFormSets:
        for descriptorValue in formset.save(commit=False):
          descriptorValue.reaction=rxn.reaction
      return redirect('reactionlist')
  else:
    reactionForm = PerformedRxnForm(user)
    reactantsFormSet = modelformset_factory(CompoundQuantity, fields=('compound', 'role', 'amount'))
    descriptorFormsets = (
      NumRxndescValFormSet(), OrdRxndescValFormSet(), BoolRxnDescValFormSet(), CatRxnDescValFormSet()
    )
  return render(request, 'reaction_form.html', {'reaction_form':reactionForm, 'reactants_formset':reactantsFormSet, 'descriptor_formsets':descriptorFormSets}) 
