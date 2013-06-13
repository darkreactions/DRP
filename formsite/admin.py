from models import Lab_Member, Lab_Group, Data
from django.contrib import admin

#class Poll_Lab_Group(admin.ModelAdmin):
	

admin.site.register(Lab_Member)
admin.site.register(Lab_Group)
admin.site.register(Data)
