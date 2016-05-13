'''A module containing views pertinent to the manipulation of reaction objects'''

import urllib
from django.views.generic import CreateView, ListView, UpdateView
from DRP.models import PerformedReaction, OrdRxnDescriptorValue, CompoundQuantity
from DRP.models import NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue
from DRP.forms import PerformedRxnForm, PerformedRxnDeleteForm
from DRP.forms import NumRxnDescValForm, OrdRxnDescValForm, BoolRxnDescValForm, CatRxnDescValForm  
from django.utils.decorators import method_decorator
from decorators import userHasLabGroup, hasSignedLicense, labGroupSelected, reactionExists
from django.contrib.auth.decorators import login_required
from django.forms.models import modelformset_factory
from DRP.forms import ModelFormSet, FormSet
from django.forms.formsets import TOTAL_FORM_COUNT
from django.shortcuts import render, redirect
from django.http import HttpResponse, Http404, HttpResponseForbidden
from django.views.decorators.http import require_POST
from django.core.exceptions import PermissionDenied
from django.contrib import messages

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
            messages.success(request, "Reaction Created Successfully")
            return redirect('editReaction', rxn.id);
    else:
        perfRxnForm = PerformedRxnForm(request.user)
    return render(request, 'reaction_create.html', {'reaction_form':perfRxnForm}) 

@login_required
@hasSignedLicense
@userHasLabGroup
@reactionExists
def addCompoundDetails(request, rxn_id):
    '''A view for adding compound details to a reaction''' 
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
            messages.success(request, 'Compound details successfully updated')
            return redirect('editReaction', rxn_id)
    else:
        formset = CompoundQuantityFormset(queryset=compoundQuantities)
    return render(request, 'reaction_detail_add.html', {'formset':formset, 'rxn_id':rxn_id, 'info_header':'Compound Quantities'})

@login_required
@hasSignedLicense
@userHasLabGroup
@reactionExists
def createGenDescVal(request, rxn_id, descValClass, descValFormClass, infoHeader):
    '''A generic view function to create descriptor values for reactions'''
    descVals = descValClass.objects.filter(reaction__id=rxn_id).filter(descriptor__calculatorSoftware="manual")
    descValFormset = modelformset_factory(model=descValClass, form=descValFormClass, can_delete=True) 
    if request.method=="POST":
        formset = descValFormset(queryset=descVals, data=request.POST)
        if formset.is_valid():
            descVals = formset.save(commit=False)
            for dv in descVals:
                dv.reaction=PerformedReactions.objects.get(id=rxn_id)
                dv.save()
            for dv in formset.deleted_objects:
                descValClass.objects.filter(id=dv.id).delete()
            messages.success(request, 'Reaction descriptor details successfully updated')
            return redirect('editReaction', rxn_id)
    else:
        formset = descValFormset(queryset=descVals)
    return render(request, 'reaction_detail_add.html', {'formset':formset, 'rxn_id':rxn_id, 'info_header':infoHeader}) 

def createTemplateRxn(request):
    pass

@login_required
@hasSignedLicense
@userHasLabGroup
@reactionExists
def editReaction(request, rxn_id):
    '''A view designed to edit performed reaction instances'''
    reaction = PerformedReaction.objects.get(id=rxn_id)
    if request.method=="POST":
        perfRxnForm = PerformedRxnForm(request.user, data=request.POST, instance=reaction)
        if perfRxnForm.is_valid():
            perfRxnForm.save()
    else:
        perfRxnForm = PerformedRxnForm(request.user, instance=reaction)
    compoundQuantities=CompoundQuantity.objects.filter(reaction__id=rxn_id)
    CompoundQuantityFormset = modelformset_factory(model=CompoundQuantity, fields=("compound", "role", "amount"), can_delete=True, extra=1)
    cqFormset = CompoundQuantityFormset(queryset=compoundQuantities) 
    descriptorFormsets = {}
    descriptorClasses=(('Numeric Descriptors', 'createNumDescVals', NumRxnDescriptorValue, NumRxnDescValForm),
                       ('Ordinal Descriptors', 'createOrdDescVals', OrdRxnDescriptorValue, OrdRxnDescValForm),
                       ('Boolean Descriptors', 'createBoolDescVals', BoolRxnDescriptorValue, BoolRxnDescValForm),
                       ('Categorical Descriptors', 'createCatDescVals', CatRxnDescriptorValue, CatRxnDescValForm))

    for descLabel, urlName, descValClass, descValFormClass in descriptorClasses: 
        descriptors = descValClass.descriptorClass.objects.filter(calculatorSoftware='manual')
        if descriptors.exists():
            descVals = descValClass.objects.filter(reaction__id=rxn_id).filter(descriptor__calculatorSoftware="manual")
            descValFormset = modelformset_factory(model=descValClass, form=descValFormClass, can_delete=True) 
            descriptorFormsets[descLabel] = {'url':urlName, 'formset': descValFormset(queryset=descVals, initial=[{'descriptor':descriptor.id} for descriptor in descriptors])}
    return render(request, 'reaction_edit.html', {'reaction_form': perfRxnForm, 'reactants_formset':cqFormset, 'descriptor_formsets':descriptorFormsets, 'reaction':reaction})

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
