"""The url configuration for the DRP project."""
from django.conf.urls import include, url
from django.contrib import admin
from . import public
from . import authentication
from . import database
from . import dashboard

admin.autodiscover()

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^', include(public.urls)),
    url(r'^', include(authentication.urls)),
    url(r'^', include(dashboard.urls)),
    url(r'^database', include(database.urls)),
]
"""The base urlconf, which includes modularised urls in the system."""
