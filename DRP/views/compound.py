"""A module containing views pertinent to compound objects."""

from django.contrib.auth.models import User
<<<<<<< HEAD
from django.views.generic import CreateView, ListView, UpdateView
from DRP.models import Compound
=======
from django.views.generic import CreateView, ListView, UpdateView, DetailView
from DRP.models import *
>>>>>>> e08a9d8bcd64b253b8f31062a7cf280d17bb3a0e
from DRP.forms import CompoundForm, LabGroupSelectionForm, CompoundEditForm, CompoundDeleteForm
from django.utils.decorators import method_decorator
from .decorators import userHasLabGroup, hasSignedLicense, labGroupSelected
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse_lazy as reverse
from django.shortcuts import render, redirect
from django.utils.http import urlencode
from django.http import HttpResponse, Http404, HttpResponseForbidden
from django.template.loader import get_template
from django.views.decorators.http import require_POST
from django.template import RequestContext


class FormFailureMixin:
    """This class adds a desireable status code to generic form views."""

    def form_invalid(self, *args, **kwargs):
        """Add a 422 status code if the form does not validate."""
        resp = super(FormFailureMixin, self).form_invalid(*args, **kwargs)
        resp.status_code = 422
        return resp


class CreateCompound(FormFailureMixin, CreateView):
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
        """Add a page heading to the normal context data."""
        context = super(CreateCompound, self).get_context_data(**kwargs)
        context['page_heading'] = 'Add a New Compound'
        return context


class EditCompound(FormFailureMixin, UpdateView):
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
            compound = Compound.objects.get(
                pk=self.get_object().pk, labGroups__in=request.user.labgroup_set.all())
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
        return redirect('compoundguide', '/')
    else:
        return HttpResponse(status=422)


class ListCompound(ListView):
    """A view managing the viewing of the compound guide."""

    template_name = 'compound_list.html'
    context_object_name = 'cg_entries'
    model = Compound

    @method_decorator(login_required)
    @method_decorator(hasSignedLicense)
    @method_decorator(userHasLabGroup)
    @labGroupSelected  # sets self.labGroup
    def dispatch(self, request, *args, **kwargs):
        """
        Overriden with a decorator to ensure that user is logged in and has at least one labGroup.

        Relate the queryset of this view to the logged in user.
        """
        self.queryset = self.labGroup.compoundguideentry_set.all()

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
            response[
                'Content-Disposition'] = 'attachment; filename="compounds.arff"'
            if 'expanded' in request.GET:
                self.queryset.toArff(response, True)
            else:
                self.queryset.toArff(response)
        else:
            raise RuntimeError(
                'The user should not be able to provoke this code')
        return response

    def get_context_data(self, **kwargs):
        """Add a lab Form and the filter formset to the existing context data."""
        context = super(ListCompound, self).get_context_data(**kwargs)
        context['lab_form'] = self.labForm
        return context
<<<<<<< HEAD
=======


@login_required
@hasSignedLicense
@userHasLabGroup
def displayCompound(request, pk):
    """Display a given compound."""
    compound = Compound.objects.get(id=pk)
    return render(request, 'compound_display.html', {'compound': compound})
>>>>>>> e08a9d8bcd64b253b8f31062a7cf280d17bb3a0e
