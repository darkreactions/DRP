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
			#Change Page
    (r'^data_transmit/(?P<num>\d+)/$', data_transmit), #Used for changing pages.
			#Upload/Download Data
    (r'^upload_CSV/$', upload_CSV),
	(r'^download_CSV/$', download_CSV),
			#Modify Data
    (r'^get_full_datum/$', get_full_datum), #"Expand" the data.
    (r'^data_update/$', data_update), #Update Information
    (r'^data_form/$', data_form), #Form for adding new data.
    (r'^compound_guide_form/$', compound_guide_form), #Form for adding new CG abbreviations.
    (r'^edit_CG_entry/$', edit_CG_entry), #Form for adding new CG abbreviations.
		#Validation
    (r'^send_CG_names/$', send_CG_names), #Send the CG name_pairs for client-side validation.
		#Authentication
    (r'^user_logout/$', user_logout), #Log Out
    (r'^user_login/$', user_login), #Log In
		#Registration
    (r'^registration_prompt/$', registration_prompt), #Redirects to correct registration choice.
    (r'^lab_registration/$', lab_registration), #Create Lab ###INACTIVE
    (r'^user_registration/$', user_registration), #Create User
	#Enable the admin:
    (r'^admin/', include(admin.site.urls)),
    
    #Send the favicon.ico:
    (r'^favicon\.ico$', 'django.views.generic.simple.redirect_to', {'url': settings.STATIC_URL+'/favicon.ico'}),
)
