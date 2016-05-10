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
    def dispatch(self, request, filetype=None, *args, **kwargs):

        if self.labGroup is not None:
                self.queryset = PerformedReaction.objects.filter(reaction_ptr__in=self.labGroup.reaction_set.all()) | PerformedReaction.objects.filter(public=True)
        else:
                self.queryset = PerformedReaction.objects.filter(public=True)

        if filetype is None:
            response = super(ListPerformedReactions, self).dispatch(request, *args, **kwargs)
        elif filetype == '.html':
            if 'page' not in request.GET:
                self.paginate_by = None
            if 'reactions_only' in request.GET:
                self.template_name='reactions_divs.html'
            response = super(ListPerformedReactions, self).dispatch(request, *args, **kwargs)
        elif filetype == '.csv':
            self.paginate_by = None
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition']='attachment; filename="reactions.csv"'
            if 'expanded' in request.GET and request.user.is_authenticated() and user.is_staff:
                self.queryset.toCsv(response, True)
            else:
                self.queryset.toCsv(response)
        elif filetype == '.arff':
            self.paginate_by = None
            response = HttpResponse(content_type='text/vnd.weka.arff')
            response['Content-Disposition']='attachment; filename="reactions.arff"'
            if 'expanded' in request.GET and request.user.is_authenticated() and user.is_staff:
                self.queryset.toArff(response, True)
            else:
                self.queryset.toArff(response)
        return response

    def get_context_data(self, **kwargs):
        context = super(ListPerformedReactions, self).get_context_data(**kwargs)
        context['lab_form'] = self.labForm
#    context['filter_formset'] = self.filterFormSet
        return context

#TODO: Remove tests for old create reaction view and 
# add ones for createReaction, editReaction, addCompoundDetails
@login_required
@hasSignedLicense
@userHasLabGroup
def createReaction(request):
    '''A view designed to create performed reaction instances'''
    if request.method == "POST":
        perfRxnForm = PerformedRxnForm(request.user, data=request.POST)
        if perfRxnForm.is_valid():
            rxn = perfRxnForm.save()
            return redirect('editReaction', rxn.id)
    else:
        perfRxnForm = PerformedRxnForm(request.user)
        return render(request, 'reaction_create.html', {'reaction_form':perfRxnForm}) 

@login_required
@hasSignedLicense
@userHasLabGroup
def addCompoundDetails(request, rxn_id):
    '''A view for adding compound details to a reaction''' 
    if PerformedReaction.objects.filter(id=rxn_id).exists() and PerformedReaction.objects.filter(labGroup__in=request.user.labgroup_set.all()):
        compoundQuantities =CompoundQuantity.objects.filter(reaction__id=rxn_id)
        CompoundQuantityFormset = modelformset_factory(model=CompoundQuantity, fields=("amount", "compound", "role"), can_delete=True)
        if request.method=="POST":
            formset = CompoundQuantityFormset(queryset=compoundQuantities, data=request.POST)
            if formset.is_valid():
                compoundQuantities = formset.save(commit=False)
                for cq in compoundQuantities:
                    cq.reaction=PerformedReaction.objects.get(id=rxn_id)
                    cq.save()
                for cq in formset.deleted_objects:
                    CompoundQuantity.objects.filter(id=cq.id).delete() #copes with a bug in deletion from django
        else:
            formset = CompoundQuantityFormset(queryset=compoundQuantities)
        return render(request, 'reaction_cq_add.html', {'formset':formset, 'rxn_id':rxn_id})
    else:
        raise Http404("This reaction cannot be found")

#TODO: Create views for each object creation separately (e.g. compound quantities) The formset objects can still be put onto the reaction edit page so long as they point to the right place.
#TODO: Create a view for creating template reactions
@login_required
@hasSignedLicense
@userHasLabGroup
def editReaction(request, pk):
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
        numRxnDescriptorValues = NumRxnDescriptorValue.objects.filter(reaction=reaction, descriptor__calculatorSoftware='manual')
        ordRxnDescriptorValues = OrdRxnDescriptorValue.objects.filter(reaction=reaction, descriptor__calculatorSoftware='manual')
        boolRxnDescriptorValues = BoolRxnDescriptorValue.objects.filter(reaction=reaction, descriptor__calculatorSoftware='manual')
        catRxnDescriptorValues = CatRxnDescriptorValue.objects.filter(reaction=reaction, descriptor__calculatorSoftware='manual') 
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
