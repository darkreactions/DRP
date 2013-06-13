from django.conf.urls.defaults import *
from formsite.views import *

# Uncomment the next two lines to enable the admin: ###C
from django.contrib import admin
admin.autodiscover()

handler500 = 'djangotoolbox.errorviews.server_error'

urlpatterns = patterns('',
    ('^_ah/warmup$', 'djangoappengine.views.warmup'),
	#Example URL:
    ###('^$', 'django.views.generic.simple.direct_to_template',{'template':'home.html'}),
    url(r'^$', data_view),
    url(r'^data_view/$', data_view), #Default Page View
    url(r'^data_view/(?P<num>\d+)/$', data_view), #Page Switch
    url(r'^upload_data/$', upload_data),
    url(r'^download_CSV/$', download_CSV),
    url(r'^data_transmit/$', data_transmit), #Better Page Switch###
    url(r'^data_transmit/(?P<num>\d+)/$', data_transmit), #Better Page Switch###
    url(r'^data_update/$', data_update), #Update Information
    url(r'^user_registration/$', user_registration), #Create User
    url(r'^user_login/$', user_login), #Log In
    url(r'^user_logout/$', user_logout), #Log Out
	#Enable the admin:
    url(r'^admin/', include(admin.site.urls)),
)
