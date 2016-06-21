"""Contains admin information for Django."""

from django.contrib import admin
from models import LabGroup, License, Compound, ChemicalClass, CompoundRole
from models import PerformedReaction, OrdRxnDescriptor, NumRxnDescriptor
from models import BoolRxnDescriptor, CatRxnDescriptor
from models import CategoricalDescriptorPermittedValue
from models import CompoundGuideEntry
from forms import LabGroupForm, CompoundAdminForm, PerformedRxnAdminForm
from forms import OrdRxnDescriptorForm, NumRxnDescriptorForm
from forms import BoolRxnDescriptorForm
from forms import CatRxnDescriptorForm, CatDescPermittedValueForm


class LabGroupAdmin(admin.ModelAdmin):

    """Specifies a form to use for administrating lab groups."""

    form = LabGroupForm


def licenseSnippet(license):
    """Return a snippet of a license."""
    return license.text[:100] + '...'

licenseSnippet.short_description = 'License Snippet'


class LicenseAdmin(admin.ModelAdmin):

    """Provides a more specific display for licenses in the django admin."""

    list_display = (licenseSnippet, 'effectiveDate')


class CompoundGuideEntryAdmin(admin.ModelAdmin):

    """Provides a more specific display for compound guide entries"""

    list_display = ('abbrev', 'labGroup', 'compound')


class CompoundAdmin(admin.ModelAdmin):

    """Specifies the form and list display of a compound in django."""

    form = CompoundAdminForm
    list_display = ('name', 'CSID', 'custom')


class ChemicalClassAdmin(admin.ModelAdmin):

    """Specifies a list display for chemical classes."""

    list_display = ('label', 'description')


class CompoundRoleAdmin(admin.ModelAdmin):

    """Specifies a list display for compound roles."""

    list_display = ('label', 'description')


class PerformedRxnAdmin(admin.ModelAdmin):

    """Defines a list display and form for Performed Reactions in django."""

    list_display = ('reference', 'user', 'labGroup', 'performedDateTime')
    form = PerformedRxnAdminForm


class BoolRxnDescriptorAdmin(admin.ModelAdmin):

    """Specifies a form and list display for custom boolean rxn descriptors."""

    list_display = ('heading', 'name')
    form = BoolRxnDescriptorForm


class NumRxnDescriptorAdmin(admin.ModelAdmin):

    """Specifies a form, list display for custom numerical rxn descriptors."""

    list_display = ('heading', 'name', 'maximum', 'minimum')
    form = NumRxnDescriptorForm

    def get_form(self, request, obj=None, **kwargs):
        """Prevent maximum and minimum values from being changed."""
        if obj is not None:
            kwargs['exclude'] = ('maximum', 'minimum')
        return super(NumRxnDescriptorAdmin, self).get_form(
            request, obj, **kwargs)


class OrdRxnDescriptorAdmin(admin.ModelAdmin):

    """Specifies a form, list display for custom numerical rxn descriptors."""

    list_display = ('heading', 'name', 'maximum', 'minimum')
    form = OrdRxnDescriptorForm

    def get_form(self, request, obj=None, **kwargs):
        """Prevent maximum and minimum values from being changed."""
        if obj is not None:
            kwargs['exclude'] = ('maximum', 'minimum')
        return super(OrdRxnDescriptorAdmin, self).get_form(
            request, obj, **kwargs)


class CatRxnDescriptorAdmin(admin.ModelAdmin):

    """Specifies a form, list display for custom category rxn descriptors."""

    list_display = ('heading', 'name')
    form = CatRxnDescriptorForm


class CatDescPermValAdmin(admin.ModelAdmin):

    """Specifies a form, list display for permitted categorical values."""

    list_display = ('value', 'descriptor')
    form = CatDescPermittedValueForm

register = admin.site.register
register(LabGroup, LabGroupAdmin)
register(License, LicenseAdmin)
register(Compound, CompoundAdmin)
register(ChemicalClass, ChemicalClassAdmin)
register(CompoundRole, CompoundRoleAdmin)
register(PerformedReaction, PerformedRxnAdmin)
register(BoolRxnDescriptor, BoolRxnDescriptorAdmin)
register(NumRxnDescriptor, NumRxnDescriptorAdmin)
register(CatRxnDescriptor, CatRxnDescriptorAdmin)
register(OrdRxnDescriptor, OrdRxnDescriptorAdmin)
register(CategoricalDescriptorPermittedValue, CatDescPermValAdmin)
register(CompoundGuideEntry, CompoundGuideEntryAdmin)
