'''The url configuration for the DRP project'''
from django.conf.urls import patterns, include, url
from django.contrib import admin
import public
import authentication
import database
import dashboard

admin.autodiscover()

urlpatterns = patterns('',
  (r'^admin/', include(admin.site.urls)),
  (r'^', include(public.urls)),
  (r'^', include(authentication.urls)),
  (r'^database/', include(database.urls)),
  #(r'^dashboard/', include(dashboard.urls))
)
'''The base urlconf, which includes modularised urls in the system'''
