from django.conf.urls import *
from django.conf import settings
from django.contrib import admin
from DRP.views import *

# Uncomment the next two lines to enable the admin: ###C
admin.autodiscover()
handler500 = 'DRP.views.display_500_error'
handler404 = 'DRP.views.display_404_error'

urlpatterns = patterns('',
  #Home and info pages.
    (r'^$', info_page, {"page":"home"}),
    (r'^home/?$', info_page, {"page":"home"}),
    (r'^papers/?$', info_page, {"page":"papers"}),
    (r'^about/?$', info_page, {"page":"about"}),
  #Database
    (r'^database/$', database), #Encompassing data view.
   #Change Page
    (r'^data_transmit/$', data_transmit), #Used for changing pages.
   #Upload/Download Data
    (r'^upload_CSV/$', upload_CSV),
    (r'^download_CSV/$', download_CSV),
   #Modify Data
    (r'^change_Data/$', change_Data), #[JSON] Edit Data Entry
    (r'^delete_Data/$', delete_Data), #[JSON] Delete Data Entries
    (r'^data_form/$', data_form), #Form for adding new data.
    (r'^compound_guide_form/$', compound_guide_form), #Form for adding new CG abbreviations.
    (r'^compound_guide_entry/$', compound_guide_entry), #Return a single CG table entry.
    (r'^edit_CG_entry/$', edit_CG_entry), #Edit a CG entry.
   #Validation
    (r'^send_CG_names/$', send_CG_names), #Send the CG name_pairs for client-side validation.
  #Predictions
    (r'^predictions/$', predictions),
    (r'^gather_SVG/$', gather_SVG),
  #Recommendations
    (r'^recommend/$', recommend),
    (r'^saved/$', saved),
    (r'^assign_user/$', assign_user_to_rec),
    (r'^save_recommendation/$', edit_recommendation, {"action":"save"}),
    (r'^unsave_recommendation/$', edit_recommendation, {"action":"unsave"}),
    (r'^sensical_recommendation/$', edit_recommendation, {"action":"sense"}),
    (r'^nonsensical_recommendation/$', edit_recommendation, {"action":"nonsense"}),
  #Searching
    (r'^search/$', search), 
  #Users and Labs
   #Authentication
    (r'^user_logout/$', user_logout), #Log Out
    (r'^user_login/$', user_login), #Log In
   #Registration
    (r'^user_license_agreement/$', get_user_license_agreement),
    (r'^update_user_license_agreement/$', update_user_license_agreement),
    (r'^registration_prompt/$', registration_prompt), #Redirects to correct registration choice.
    (r'^lab_registration/$', lab_registration), #Create Lab ###INACTIVE
    (r'^user_registration/$', user_registration), #Create User
    (r'^user_update/$', user_update), #Create User
    (r'^change_password/$', change_password), #Change Password

 #Enable the admin:
    (r'^admin/', include(admin.site.urls)),
)
