from django.contrib import admin
from models import LabGroup, License
from forms import LabGroupForm

class LabGroupAdmin(admin.ModelAdmin):
  form = LabGroupForm


admin.site.register(LabGroup, LabGroupAdmin)
admin.site.register(License)
