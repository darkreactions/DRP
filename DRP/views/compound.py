"""A module containing views pertinent to compound objects."""

from django.contrib.auth.models import User
from django.views.generic import CreateView, ListView, UpdateView
from DRP.models import Compound
from DRP.forms import CompoundForm, LabGroupSelectionForm, CompoundEditForm, CompoundDeleteForm, CompoundUploadForm, CompoundFilterForm
from DRP.forms import CompoundFilterFormSet, AdvancedCompoundFilterFormSet
from django.utils.decorators import method_decorator
from decorators import userHasLabGroup, hasSignedLicense, labGroupSelected
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse_lazy as reverse
from django.shortcuts import render, redirect
from django.utils.http import urlencode
from django.http import HttpResponse, Http404, HttpResponseForbidden
from django.template.loader import get_template
from django.views.decorators.http import require_POST
from django.template import RequestContext


class CreateCompound(CreateView):

    """A view managing the creation of compound objects."""

    model = Compound
    form_class = CompoundForm
    template_name = 'compound_form.html'
    success_url = reverse('compoundguide', args=['/'])

    def get_form_kwargs(self):
        """Overridden to add the request.user value into the kwargs."""
        kwargs = super(CreateCompound, self).get_form_kwargs()
        kwargs['user'] = self.request.user
        return kwargs

    @method_decorator(login_required)
    @method_decorator(hasSignedLicense)
    @method_decorator(userHasLabGroup)
    def dispatch(self, request, *args, **kwargs):
        """Overridden with a decorator to ensure that a user is at least logged in."""
        return super(CreateCompound, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(CreateCompound, self).get_context_data(**kwargs)
        context['page_heading'] = 'Add a New Compound'
        return context


class EditCompound(UpdateView):

    """A view managing the editing of compound objects."""

    form_class = CompoundEditForm
    template_name = 'compound_edit.html'
    success_url = reverse('compoundguide', args=['/'])
    model = Compound

    @method_decorator(login_required)
    @method_decorator(hasSignedLicense)
    @method_decorator(userHasLabGroup)
    def dispatch(self, request, *args, **kwargs):
        """Check user has sufficient credentials and has row-level permissions for this compound."""
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
    """A view managing the deletion of compound objects."""
    form = CompoundDeleteForm(data=request.POST, user=request.user)
    if form.is_valid():
        form.save()
    else:
        raise RuntimeError(str(form.errors))
    return redirect('compoundguide', '/')


@login_required
@hasSignedLicense
@userHasLabGroup
def uploadCompound(request, *args, **kwargs):
    """A view managing the upload of compound csvs."""
    if request.method == 'POST':
        form = CompoundUploadForm(data=request.POST, files=request.FILES, user=request.user)
        if form.is_valid():
            form.save()
            return redirect('compoundguide', '/')
        else:
            return render(request, 'compound_upload.html', {'form': form})
    else:
        form = CompoundUploadForm(user=request.user)
        return render(request, 'compound_upload.html', {'form': form})


class ListCompound(ListView):

    """A view managing the viewing of the compound guide."""

    template_name = 'compound_list.html'
    context_object_name = 'compounds'
    model = Compound
    formSetClass = CompoundFilterFormSet

    @method_decorator(login_required)
    @method_decorator(hasSignedLicense)
    @method_decorator(userHasLabGroup)
    @labGroupSelected  # sets self.labGroup
    def dispatch(self, request, *args, **kwargs):
        """
        Overriden with a decorator to ensure that user is logged in and has at least one labGroup.

        Relate the queryset of this view to the logged in user.
        """

        self.queryset = self.labGroup.compound_set.all()

        if 'filter' in request.GET:
            self.filterFormSet = self.formSetClass(user=request.user, labGroup=self.labGroup, data=request.GET)
            if self.filterFormSet.is_valid():
                self.queryset = self.filterFormSet.fetch()
                self.filterFormSet = self.formSetClass(user=request.user, labGroup=self.labGroup, initial=self.filterFormSet.cleaned_data)
            else:
                self.queryset = Compound.objects.none()
        else:
            self.filterFormSet = self.formSetClass(user=request.user, labGroup=self.labGroup)

        fileType = kwargs.get('filetype')

        if fileType in ('/', '.html', None):
            return super(ListCompound, self).dispatch(request, *args, **kwargs)
        elif fileType == '.csv':
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = 'attachment; filename="compounds.csv"'
            if 'expanded' in request.GET:
                self.queryset.toCsv(response, True)
            else:
                self.queryset.toCsv(response)
        elif fileType == '.arff':
            response = HttpResponse(content_type='text/vnd.weka.arff')
            response['Content-Disposition'] = 'attachment; filename="compounds.arff"'
            if 'expanded' in request.GET:
                self.queryset.toArff(response, True)
            else:
                self.queryset.toArff(response)
        else:
            raise RuntimeError('The user should not be able to provoke this code')
        return response

    def get_context_data(self, **kwargs):
        context = super(ListCompound, self).get_context_data(**kwargs)
        context['lab_form'] = self.labForm
        context['filter_formset'] = self.filterFormSet
        return context


class AdvancedCompoundSearchView(ListCompound):

    """View for the advanced compounds search."""

    formSetClass = AdvancedCompoundFilterFormSet
    template_name = 'adv_compound_search.html'
