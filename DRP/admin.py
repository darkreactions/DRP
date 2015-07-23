from django.contrib import admin
from models.LabGroup import LabGroup
from forms import LabGroupForm

class LabGroupAdmin(admin.ModelAdmin):
  form = LabGroupForm

admin.site.register(LabGroup, LabGroupAdmin)
