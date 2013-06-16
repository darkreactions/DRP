from django.conf.urls.defaults import *
from formsite.views import *
from django.conf import settings
from django.conf.urls.static import static

# Uncomment the next two lines to enable the admin: ###C
from django.contrib import admin
admin.autodiscover()



handler500 = 'formsite.views.display_500_error'
handler404 = 'formsite.views.display_404_error'

urlpatterns = patterns('',
    ('^_ah/warmup$', 'djangoappengine.views.warmup'),
	#Example URL:
    ###('^$', 'django.views.generic.simple.direct_to_template',{'template':'home.html'}),
    (r'^$', data_view),
    (r'^data_view/$', data_view), #Default Page View
    (r'^data_view/(?P<num>\d+)/$', data_view), #Page Switch
    (r'^upload_data/$', upload_data),
    (r'^download_CSV/$', download_CSV),
    (r'^data_transmit/$', data_transmit), #Better Page Switch###
    (r'^data_transmit/(?P<num>\d+)/$', data_transmit), #Better Page Switch###
    (r'^data_update/$', data_update), #Update Information
    (r'^user_registration/$', user_registration), #Create User
    (r'^user_login/$', user_login), #Log In
    (r'^user_logout/$', user_logout), #Log Out
	#Enable the admin:
    (r'^admin/', include(admin.site.urls)),
    
    #Send the favicon.ico:
    (r'^favicon\.ico$', 'django.views.generic.simple.redirect_to', {'url': settings.STATIC_URL+'favicon.ico'}),
)
