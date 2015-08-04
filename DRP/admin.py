from django.contrib import admin
from models import LabGroup, License
from forms import LabGroupForm

class LabGroupAdmin(admin.ModelAdmin):
  form = LabGroupForm

def licenseSnippet(license):
  return license.text[:100] + '...'

licenseSnippet.short_description = 'License Snippet'

class LicenseAdmin(admin.ModelAdmin):

  list_display = (licenseSnippet, 'effectiveDate')

admin.site.register(LabGroup, LabGroupAdmin)
admin.site.register(License, LicenseAdmin)
