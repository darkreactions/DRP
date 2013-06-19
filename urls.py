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
		#Database
    (r'^$', database),
    (r'^database/$', database), #Encompassing data view.
    (r'^data_form/$', data_form), #Encompassing data view.
    (r'^upload_CSV/$', upload_CSV),
		#Authentication
    (r'^user_logout/$', user_logout), #Log Out
    (r'^user_login/$', user_login), #Log In
    (r'^registration_prompt/$', registration_prompt), #Redirects to correct registration choice.
    (r'^lab_registration/$', lab_registration), #Create Lab ###INACTIVE
    (r'^user_registration/$', user_registration), #Create User
    #
    (r'^download_CSV/$', download_CSV),
    (r'^data_transmit/$', data_transmit), #Better Page Switch###
    (r'^data_transmit/(?P<num>\d+)/$', data_transmit), #Better Page Switch###
    (r'^data_update/$', data_update), #Update Information
	#Enable the admin:
    (r'^admin/', include(admin.site.urls)),
    
    #Send the favicon.ico:
    (r'^favicon\.ico$', 'django.views.generic.simple.redirect_to', {'url': settings.STATIC_URL+'/favicon.ico'}),
)
