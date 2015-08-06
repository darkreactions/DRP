from django.contrib import admin
from models import LabGroup, License, Compound, ChemicalClass 
from forms import LabGroupForm, CompoundAdminForm

class LabGroupAdmin(admin.ModelAdmin):
  form = LabGroupForm

def licenseSnippet(license):
  return license.text[:100] + '...'

licenseSnippet.short_description = 'License Snippet'

class LicenseAdmin(admin.ModelAdmin):

  list_display = (licenseSnippet, 'effectiveDate')

class CompoundAdmin(admin.ModelAdmin):
  form = CompoundAdminForm

admin.site.register(LabGroup, LabGroupAdmin)
admin.site.register(License, LicenseAdmin)
admin.site.register(Compound, CompoundAdmin)
admin.site.register(ChemicalClass)
