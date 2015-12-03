'''A module containing views pertinent to the manipulation of reaction objects'''

from django.views.generic import CreateView, ListView, UpdateView
from DRP.models import PerformedReaction, OrdRxnDescriptorValue, CompoundQuantity
from DRP.models import NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue
from DRP.forms import PerformedRxnForm, PerformedRxnDeleteForm
from DRP.forms import NumRxnDescValForm, OrdRxnDescValForm, BoolRxnDescValForm, CatRxnDescValForm  
from django.utils.decorators import method_decorator
from decorators import userHasLabGroup, hasSignedLicense, labGroupSelected
from django.contrib.auth.decorators import login_required
from django.forms.models import modelformset_factory
from DRP.forms import ModelFormSet, FormSet
from django.forms.formsets import TOTAL_FORM_COUNT
from django.shortcuts import render, redirect
from django.http import HttpResponse, Http404, HttpResponseForbidden
from django.views.decorators.http import require_POST
from django.core.exceptions import PermissionDenied

class ListPerformedReactions(ListView):
  '''Standard list view of performed reactions, adjusted to deal with a few DRP idiosyncrasies'''

  template_name='reactions_list.html'
  context_object_name='reactions'
  model=PerformedReaction 
  paginate_by=20
  
  @labGroupSelected #sets self.labGroup
  def dispatch(self, request, *args, **kwargs):

    if self.labGroup is not None:
        self.queryset = PerformedReaction.objects.filter(reaction_ptr__in=self.labGroup.reaction_set.all()) | PerformedReaction.objects.filter(public=True)
    else:
        self.queryset = PerformedReaction.objects.filter(public=True)

    fileType = kwargs.get('filetype')
    if fileType in ('/', '.html', None):
      return super(ListPerformedReactions, self).dispatch(request, *args, **kwargs)
    elif fileType == '.csv':
      self.paginate_by = None
      response = HttpResponse(content_type='text/csv')
      response['Content-Disposition']='attachment; filename="compounds.csv"'
      if 'expanded' in request.GET and request.user.is_authenticated():
        self.queryset.toCsv(response, True)
      else:
        self.queryset.toCsv(response)
    elif fileType == '.arff':
      self.paginate_by = None
      response = HttpResponse(content_type='text/vnd.weka.arff')
      response['Content-Disposition']='attachment; filename="compounds.arff"'
      if 'expanded' in request.GET and request.user.is_authenticated():
        self.queryset.toArff(response, True)
      else:
        self.queryset.toArff(response)
    else:
      raise RuntimeError('The user should not be able to provoke this code')
    return response

  def get_context_data(self, **kwargs):
    context = super(ListPerformedReactions, self).get_context_data(**kwargs)
    context['lab_form'] = self.labForm
#    context['filter_formset'] = self.filterFormSet
    return context

@login_required
@hasSignedLicense
@userHasLabGroup
def reactionForm(request, pk=None):
  '''A view designed to create performed reaction instances'''
  descFields = ('descriptor', 'value')
  if pk == None:
    reaction = None
  else:
    try:
      reaction=PerformedReaction.objects.get(pk=pk)
    except PerformedReaction.DoesNotExist:
      raise Http404("This reaction cannot be found")
  if reaction is not None:
    reactants = CompoundQuantity.objects.filter(reaction=reaction.reaction_ptr)
    numRxnDescriptorValues = NumRxnDescriptorValue.objects.filter(reaction=reaction.reaction_ptr)
    ordRxnDescriptorValues = OrdRxnDescriptorValue.objects.filter(reaction=reaction.reaction_ptr)
    boolRxnDescriptorValues = BoolRxnDescriptorValue.objects.filter(reaction=reaction.reaction_ptr)
    catRxnDescriptorValues = CatRxnDescriptorValue.objects.filter(reaction=reaction.reaction_ptr) 
  else:
    reactants=None
    numRxnDescriptorValues = None 
    ordRxnDescriptorValues = None
    boolRxnDescriptorValues = None
    catRxnDescriptorValues = None
    
  if request.method=='POST':
    reactantsFormSetInst = ModelFormSet(CompoundQuantity, fields=('compound', 'role', 'amount'), data=request.POST, canAdd=True, canDelete=True, instances=reactants)
    reactionForm = PerformedRxnForm(request.user, data=request.POST, instance=reaction) 

    descriptorFormSets = (
      ModelFormSet(NumRxnDescriptorValue, formClass=NumRxnDescValForm, data=request.POST, prefix='num', canDelete=True, instances=numRxnDescriptorValues),
      ModelFormSet(OrdRxnDescriptorValue, formClass=OrdRxnDescValForm, data=request.POST, prefix='ord', canDelete=True, instances=ordRxnDescriptorValues),
      ModelFormSet(BoolRxnDescriptorValue, formClass=BoolRxnDescValForm, data=request.POST, prefix='bool', canDelete=True, instances=boolRxnDescriptorValues),
      ModelFormSet(CatRxnDescriptorValue, formClass=CatRxnDescValForm, data=request.POST, prefix='cat', canDelete=True, instances=catRxnDescriptorValues)
    )

    if 'save' in request.POST: 
      if reactionForm.is_valid() and reactantsFormSetInst.is_valid() and all(d.is_valid() for d in descriptorFormSets):
        rxn = reactionForm.save()
        reactants = reactantsFormSetInst.save(commit=False)
        for reactant in reactants:
          reactant.reaction=rxn.reaction_ptr
          reactant.save()
        CompoundQuantity.objects.filter(reaction=rxn.reaction_ptr).exclude(pk__in=(reactant.pk for reactant in reactants)).delete()
        cdvs = [] #collated descriptor values
        for formSet in descriptorFormSets:
          descriptorValues = formSet.save(commit=False)
          cdvs.append(descriptorValues)
          for descriptorValue in descriptorValues:
            descriptorValue.reaction=rxn.reaction_ptr
            descriptorValue.save()
        NumRxnDescriptorValue.objects.filter(reaction=rxn.reaction_ptr).exclude(pk__in=(dv.pk for dv in cdvs[0]))
        OrdRxnDescriptorValue.objects.filter(reaction=rxn.reaction_ptr).exclude(pk__in=(dv.pk for dv in cdvs[1]))
        BoolRxnDescriptorValue.objects.filter(reaction=rxn.reaction_ptr).exclude(pk__in=(dv.pk for dv in cdvs[2]))
        CatRxnDescriptorValue.objects.filter(reaction=rxn.reaction_ptr).exclude(pk__in=(dv.pk for dv in cdvs[3]))
        return redirect('reactionlist')
  else:
    reactionForm = PerformedRxnForm(request.user, instance=reaction)
    reactantsFormSetInst = ModelFormSet(CompoundQuantity, fields=('compound', 'role', 'amount'), canAdd=True, instances=reactants, canDelete=True)
    descriptorFormSets = (
      ModelFormSet(NumRxnDescriptorValue, formClass=NumRxnDescValForm, prefix='num', instances=numRxnDescriptorValues, canDelete=True),
      ModelFormSet(OrdRxnDescriptorValue, formClass=OrdRxnDescValForm, prefix='ord', instances=ordRxnDescriptorValues, canDelete=True),
      ModelFormSet(BoolRxnDescriptorValue, formClass=BoolRxnDescValForm, prefix='bool', instances=boolRxnDescriptorValues, canDelete=True),
      ModelFormSet(CatRxnDescriptorValue, formClass=CatRxnDescValForm, prefix='cat', instances=catRxnDescriptorValues, canDelete=True)
    )
  return render(request, 'reaction_form.html', {'reaction_form':reactionForm, 'reactants_formset':reactantsFormSetInst, 'descriptor_formsets':descriptorFormSets}) 

@require_POST
@login_required
@hasSignedLicense
@userHasLabGroup
def deleteReaction(request, *args, **kwargs):
  '''A view managing the deletion of reaction objects'''
  form =PerformedRxnDeleteForm(data=request.POST, user=request.user) 
  if form.is_valid():
    form.save()
  else:
    raise RuntimeError(str(form.errors))
  return redirect('reactionlist')

@require_POST
@login_required
@hasSignedLicense
@userHasLabGroup
def invalidateReaction(request, *args, **kwargs):
  '''A view managing the deletion of reaction objects'''
  form =PerformedRxnInvalidateForm(data=request.POST, user=request.user) 
  if form.is_valid():
    form.save()
  else:
    raise RuntimeError(str(form.errors))
  return redirect('reactionlist')
