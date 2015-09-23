from django.contrib import admin
from models import LabGroup, License, Compound, ChemicalClass, CompoundRole
from models import PerformedReaction, OrdRxnDescriptor, NumRxnDescriptor
from models import BoolRxnDescriptor, CatRxnDescriptor, CategoricalDescriptorPermittedValue 
from forms import LabGroupForm, CompoundAdminForm, PerformedRxnAdminForm
from forms import OrdRxnDescriptorForm, NumRxnDescriptorForm, BoolRxnDescriptorForm
from forms import CatRxnDescriptorForm, CatDescPermittedValueForm

class LabGroupAdmin(admin.ModelAdmin):
    form = LabGroupForm

def licenseSnippet(license):
    return license.text[:100] + '...'

licenseSnippet.short_description = 'License Snippet'

class LicenseAdmin(admin.ModelAdmin):

    list_display = (licenseSnippet, 'effectiveDate')

class CompoundAdmin(admin.ModelAdmin):
    form = CompoundAdminForm

    list_display = ('abbrev', 'name', 'CSID', 'custom', 'labGroup')

class ChemicalClassAdmin(admin.ModelAdmin):
    
    list_display = ('label', 'description')

class CompoundRoleAdmin(admin.ModelAdmin):

    list_display = ('label', 'description')

class PerformedRxnAdmin(admin.ModelAdmin):

    list_display = ('reference', 'user', 'labGroup', 'performedDateTime')
    form = PerformedRxnAdminForm

class BoolRxnDescriptorAdmin(admin.ModelAdmin):

    list_display=('heading', 'name')
    form=BoolRxnDescriptorForm

class NumRxnDescriptorAdmin(admin.ModelAdmin):
    
    list_display=('heading', 'name', 'maximum', 'minimum')
    form=NumRxnDescriptorForm

    def get_form(self, request, obj=None, **kwargs):
        if obj is not None:
            kwargs['exclude'] = ('maximum', 'minimum')
        return super(NumRxnDescriptorAdmin, self).get_form(request, obj, **kwargs)

class OrdRxnDescriptorAdmin(admin.ModelAdmin): 

    list_display=('heading', 'name', 'maximum', 'minimum')
    form=OrdRxnDescriptorForm

    def get_form(self, request, obj=None, **kwargs):
        if obj is not None:
            kwargs['exclude'] = ('maximum', 'minimum')
        return super(OrdRxnDescriptorAdmin, self).get_form(request, obj, **kwargs)

class CatRxnDescriptorAdmin(admin.ModelAdmin):
    list_display=('heading', 'name')
    form = CatRxnDescriptorForm

class CatDescPermValAdmin(admin.ModelAdmin):
    list_display=('value', 'descriptor')
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
