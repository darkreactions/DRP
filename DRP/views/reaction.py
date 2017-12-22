"""A module containing views pertinent to the manipulation of reaction objects."""
from django.conf import settings
from django.views.generic import CreateView, ListView, UpdateView
from DRP.models import PerformedReaction, OrdRxnDescriptorValue, CompoundQuantity
from DRP.models import NumRxnDescriptorValue, BoolRxnDescriptorValue, CatRxnDescriptorValue
from DRP.forms import PerformedRxnForm, PerformedRxnDeleteForm
from DRP.forms import NumRxnDescValFormFactory, OrdRxnDescValFormFactory, BoolRxnDescValFormFactory, CatRxnDescValFormFactory
from DRP.forms import PerformedRxnInvalidateForm, PerformedRxnDeleteForm
from django.utils.decorators import method_decorator
from .decorators import userHasLabGroup, hasSignedLicense, labGroupSelected, reactionExists
from django.contrib.auth.decorators import login_required
from django.forms.models import modelformset_factory
from DRP.forms import compoundQuantityFormFactory
from django.forms.formsets import TOTAL_FORM_COUNT
from django.shortcuts import render
from .helpers import redirect
from django.http import HttpResponse, Http404, HttpResponseForbidden
from django.views.decorators.http import require_POST
from django.core.exceptions import PermissionDenied
from django.contrib import messages
import logging


@login_required
@hasSignedLicense
@userHasLabGroup
@reactionExists
def labBookImage(request, labgroup_id, reference=None):
    """view that does security checking before prompting the server to show an image."""
    response = HttpResponse()
    response[settings.MEDIA_X_HEADER] = settings.SECURE_MEDIA_URL + '/lab_notes/' + \
        PerformedReaction.objects.get(
            reference=reference, labGroup__id=labgroup_id).labBookPage.name
    response['Content-Type'] = 'image/jpeg'
    return response


class ListPerformedReactions(ListView):
    """Standard list view of performed reactions, adjusted to deal with a few DRP idiosyncrasies."""

    template_name = 'reactions_list.html'
    context_object_name = 'reactions'
    model = PerformedReaction
    paginate_by = 20

    @labGroupSelected  # sets self.labGroup
    def dispatch(self, request, filetype=None, *args, **kwargs):
        """Render the view according to the appropriate filetype."""
        if self.labGroup is not None:
            self.queryset = PerformedReaction.objects.filter(
                reaction_ptr__in=self.labGroup.reaction_set.all()) | PerformedReaction.objects.filter(public=True)
        else:
            self.queryset = PerformedReaction.objects.filter(public=True)
        self.queryset = self.queryset.order_by('-insertedDateTime')
        if filetype is None:
            response = super(ListPerformedReactions, self).dispatch(
                request, *args, **kwargs)
        elif filetype == '.html':
            if 'page' not in request.GET:
                self.paginate_by = None
            if 'reactions_only' in request.GET:
                self.template_name = 'reactions_divs.html'
            response = super(ListPerformedReactions, self).dispatch(
                request, *args, **kwargs)
        elif filetype == '.csv':
            self.paginate_by = None
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = 'attachment; filename="reactions.csv"'
            if 'expanded' in request.GET and request.user.is_authenticated() and request.user.is_staff:
                headers = self.load_from_dsc(settings.RECOMMENDED_WHITE_LIST)
                self.queryset.toCsv(response, expanded=True, whitelistHeaders=headers)
            else:
                self.queryset.toCsv(response)
        elif filetype == '.arff':
            self.paginate_by = None
            response = HttpResponse(content_type='text/vnd.weka.arff')
            response[
                'Content-Disposition'] = 'attachment; filename="reactions.arff"'
            if 'expanded' in request.GET and request.user.is_authenticated() and request.user.is_staff:
                self.queryset.toArff(response, True)
            else:
                self.queryset.toArff(response)
        return response

    def get_context_data(self, **kwargs):
        """Attach the lab form as additional context; deprecated."""
        context = super(ListPerformedReactions,
                        self).get_context_data(**kwargs)
        context['lab_form'] = self.labForm
        return context

    def load_from_dsc(self, dsc_file):
        """Load from dsc file."""
        with open(dsc_file, 'r') as f:
            headers = f.read().split('\n')
            headers = headers[:len(headers) - 1]
        return headers


@login_required
@hasSignedLicense
@userHasLabGroup
def createReaction(request):
    """A view designed to create performed reaction instances."""
    status = 200
    if request.method == "POST":
        perfRxnForm = PerformedRxnForm(
            request.user, data=request.POST, files=request.FILES)
        if perfRxnForm.is_valid():
            rxn = perfRxnForm.save()
            messages.success(request, "Reaction Created Successfully")
            return redirect('addCompoundDetails', rxn.id, params={'creating': True})
        else:
            status = 422
    else:
        perfRxnForm = PerformedRxnForm(request.user)
    return render(request, 'reaction_create.html', {'reaction_form': perfRxnForm}, status=status)


@login_required
@hasSignedLicense
@userHasLabGroup
@reactionExists
def addCompoundDetails(request, rxn_id):
    """A view for adding compound details to a reaction."""
    status = 200
    compoundQuantities = CompoundQuantity.objects.filter(reaction__id=rxn_id)
    CompoundQuantityFormset = modelformset_factory(model=CompoundQuantity, form=compoundQuantityFormFactory(
        rxn_id), can_delete=('creating' not in request.GET), extra=6)
    if request.method == "POST":
        formset = CompoundQuantityFormset(
            queryset=compoundQuantities, data=request.POST, prefix='quantities')
        if formset.is_valid():
            formset.save()
            # copes with a bug in deletion from django
            CompoundQuantity.objects.filter(
                id__in=[cq.id for cq in formset.deleted_objects]).delete()
            messages.success(request, 'Compound details successfully updated')
            if 'creating' in request.GET:
                return redirect('createNumDescVals', rxn_id, params={'creating': True})
            else:
                return redirect('editReaction', rxn_id)
        else:
            status = 422
    else:
        formset = CompoundQuantityFormset(
            queryset=compoundQuantities, prefix='quantities')
    return render(request, 'reaction_compound_add.html', {'reactants_formset': formset, 'reaction': PerformedReaction.objects.get(id=rxn_id), }, status=status)


@login_required
@hasSignedLicense
@userHasLabGroup
@reactionExists
def createGenDescVal(request, rxn_id, descValClass, descValFormClass, infoHeader, createNext):
    """A generic view function to create descriptor values for reactions."""
    descVals = descValClass.objects.filter(reaction__id=rxn_id).filter(
        descriptor__calculatorSoftware="manual")
    descriptors = descValClass.descriptorClass.objects.filter(
        calculatorSoftware="manual")
    initialDescriptors = descValClass.descriptorClass.objects.filter(
        isDefaultForLabGroups__reaction__id=rxn_id,
        calculatorSoftware='manual').exclude(id__in=set(descVal.descriptor.id for descVal in descVals))
    descValFormset = modelformset_factory(model=descValClass, form=descValFormClass(
        rxn_id), can_delete=('creating' not in request.GET), extra=initialDescriptors.count())
    status = 200
    if ('creating' in request.GET and initialDescriptors.exists()) or descriptors.exists():
        if request.method == "POST":
            formset = descValFormset(
                queryset=descVals, data=request.POST, prefix=request.resolver_match.url_name)
            # this weird prefix is caused by the generic nature of this function and the neccessity to namespace
            # the different form elements in the formsets used in the edit
            # reaction view.
            if formset.is_valid():
                descVals = formset.save()
                descValClass.objects.filter(
                    pk__in=[dv.pk for dv in formset.deleted_objects]).delete()
                messages.success(
                    request, 'Reaction descriptor details successfully updated')
                if createNext is None or 'creating' not in request.GET:
                    return redirect('editReaction', rxn_id)
                else:
                    return redirect(createNext, rxn_id, params={'creating': True})
            else:
                status = 200
        else:
            formset = descValFormset(queryset=descVals, initial=[
                                     {'descriptor': descriptor.id} for descriptor in initialDescriptors], prefix=request.resolver_match.url_name)
        return render(request, 'reaction_detail_add.html', {'formset': formset, 'rxn_id': rxn_id, 'info_header': infoHeader}, status=status)
    elif createNext is None or 'creating' not in request.GET:
        return redirect('editReaction', rxn_id)
    else:
        return redirect(createNext, rxn_id, params={'creating': True})


@login_required
@hasSignedLicense
@userHasLabGroup
@reactionExists
def editReaction(request, rxn_id):
    """A view designed to edit performed reaction instances."""
    status = 200
    reaction = PerformedReaction.objects.get(id=rxn_id)
    if request.method == "POST":
        perfRxnForm = PerformedRxnForm(
            request.user, data=request.POST, instance=reaction, files=request.FILES)
        if perfRxnForm.is_valid():
            perfRxnForm.save()
            messages.success(request, "Reaction successfully updated.")
        else:
            status = 422
    else:
        perfRxnForm = PerformedRxnForm(request.user, instance=reaction)
    compoundQuantities = CompoundQuantity.objects.filter(reaction__id=rxn_id)
    CompoundQuantityFormset = modelformset_factory(
        model=CompoundQuantity, form=compoundQuantityFormFactory(rxn_id), can_delete=True, extra=1)
    cqFormset = CompoundQuantityFormset(
        queryset=compoundQuantities, prefix="quantities")
    descriptorFormsets = {}
    descriptorClasses = (('Numeric Descriptors', 'createNumDescVals', NumRxnDescriptorValue, NumRxnDescValFormFactory(rxn_id)),
                         ('Ordinal Descriptors', 'createOrdDescVals',
                          OrdRxnDescriptorValue, OrdRxnDescValFormFactory(rxn_id)),
                         ('Boolean Descriptors', 'createBoolDescVals',
                          BoolRxnDescriptorValue, BoolRxnDescValFormFactory(rxn_id)),
                         ('Categorical Descriptors', 'createCatDescVals', CatRxnDescriptorValue, CatRxnDescValFormFactory(rxn_id)))
    for descLabel, urlName, descValClass, descValFormClass in descriptorClasses:
        descVals = descValClass.objects.filter(reaction__id=rxn_id).filter(
            descriptor__calculatorSoftware="manual")
        descriptors = descValClass.descriptorClass.objects.filter(
            calculatorSoftware='manual')
        initialDescriptors = descValClass.descriptorClass.objects.filter(
            isDefaultForLabGroups=reaction.labGroup,
            calculatorSoftware='manual').exclude(id__in=set(descVal.descriptor.id for descVal in descVals))
        if descriptors.exists():
            descValFormset = modelformset_factory(
                model=descValClass, form=descValFormClass, can_delete=True, extra=initialDescriptors.count() + 1)
            descriptorFormsets[descLabel] = {'url': urlName, 'formset': descValFormset(queryset=descVals, initial=[
                                                                                       {'descriptor': descriptor.id} for descriptor in initialDescriptors], prefix=urlName)}
    return render(request, 'reaction_edit.html', {'reaction_form': perfRxnForm, 'reactants_formset': cqFormset, 'descriptor_formsets': descriptorFormsets, 'reaction': reaction}, status=status)


@require_POST
@login_required
@hasSignedLicense
@userHasLabGroup
def deleteReaction(request, *args, **kwargs):
    """A view managing the deletion of reaction objects."""
    form = PerformedRxnDeleteForm(data=request.POST, user=request.user)
    if form.is_valid():
        form.save()
        return redirect('reactionlist')
    else:
        return HttpResponse(status=422)


@require_POST
@login_required
@hasSignedLicense
@userHasLabGroup
def invalidateReaction(request, *args, **kwargs):
    """A view managing the deletion of reaction objects."""
    form = PerformedRxnInvalidateForm(data=request.POST, user=request.user)
    if form.is_valid():
        form.save()
        return redirect('reactionlist')
    else:
        return HttpResponse(status=422)
