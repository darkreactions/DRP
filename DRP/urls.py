from django.conf.urls import *
from django.conf import settings
from django.contrib import admin
from core_views import *

# Uncomment the next two lines to enable the admin: ###C
admin.autodiscover()
handler500 = 'DRP.views.errors.display_500_error'
handler404 = 'DRP.views.errors.display_404_error'

urlpatterns = patterns('',
 #Home and info pages.
    (r'^$', info_page, {"page":"home"}),
    (r'^home/?$', info_page, {"page":"home"}),
    (r'^papers/?$', info_page, {"page":"papers"}),
    (r'^about/?$', global_page, {"page": "about"}),
    (r'^contact/?$', info_page, {"page":"contact"}),
    (r'^contact_form/?$', "DRP.views.contact.contact_form"),
    (r'^explore/?$', global_page, {"page":"explore"}),
    #(r'^help/?$', info_page, {"page":"help"}),
  #Dashboard
    url(r'^dashboard/?$', "DRP.views.dashboard.get_dashboard"), #Displays the empty dashboard.
    url(r'^get_stats/?$', "DRP.views.dashboard.get_stats_json"), #Actually loads the json.
  #Database
    (r'^database/$', database), #Encompassing data view.
   #Change Page
    (r'^data_transmit/$', data_transmit), #Used for changing pages.
    (r'^recommend_transmit/$', recommendation_transmit), #Used for changing pages
   #Upload/Download database.
    (r'^upload_prompt/$', upload_prompt),
    (r'^upload_data/$', upload_CSV),
    (r'^download_prompt/$', "DRP.views.download.download_prompt"),
    (r'^download_data/$', "DRP.views.download.download_CSV"),
   #Modify Data
    (r'^change_Data/$', change_Data), #[JSON] Edit Data Entry
    (r'^delete_Data/$', delete_Data), #[JSON] Delete Data Entries
    (r'^add_reactant/$', add_reactant), #[JSON] Add Reactant Group
    (r'^delete_reactant/$', delete_reactant), #[JSON] Delete Reactant Group
    (r'^data_form/$', data_form), #Form for adding new data.
   #Add Compound Entries
    (r'^check_compound/$', check_compound), #Form for adding new CG abbreviations.
    (r'^compound_guide/$', compound_guide), #Form for adding new CG abbreviations.
    (r'^compound_guide_form/$', compound_guide_form), #Form for adding new CG abbreviations.
    (r'^compound_guide_entry/$', compound_guide_entry), #Return a single CG table entry.
    (r'^edit_CG_entry/$', edit_CG_entry), #Edit a CG entry.
   #Validation
    (r'^send_CG_names/$', "DRP.views.jsonViews.send_CG_names"), #Send the CG name_pairs for client-side validation.
  #Visualization
    (r'^get_graph/$', "DRP.views.explore_vis.get_graph_data"),
    (r'^setup_graph/$', "DRP.views.explore_vis.store_graph"),
  #Recommendations
    (r'^make_seed_recommendations/$', "DRP.views.seed_recommend.make_seed_recommendations"),
    (r'^seed/$', "DRP.views.seed_recommend.seed_recommend"),
    (r'^seed_recommend/$', "DRP.views.seed_recommend.seed_recommend"),
    (r'^check_seed_oven/$', "DRP.views.seed_recommend.check_seed_worker_cache"),
    (r'^recommend/$', recommend),
    (r'^saved/$', saved),
    (r'^change_Recommendation/$', change_Recommendation), #[JSON] Edit Rec Entry
    (r'^assign_user/$', assign_user_to_rec),
    (r'^show_recommendation/$', edit_recommendation, {"action":"show"}),
    (r'^hide_recommendation/$', edit_recommendation, {"action":"hide"}),
    (r'^save_recommendation/$', edit_recommendation, {"action":"save"}),
    (r'^unsave_recommendation/$', edit_recommendation, {"action":"unsave"}),
    (r'^sensical_recommendation/$', edit_recommendation, {"action":"sense"}),
    (r'^nonsensical_recommendation/$', edit_recommendation, {"action":"nonsense"}),
  #Rank Recommendations
    (r'^rank/$', rank),
    (r'^send_and_receive_rank/$', send_and_receive_rank),
  #Searching
    (r'^search/$', search, {"model":"Data"}),
    (r'^search/Data/?$', search, {"model":"Data"}),
    (r'^search/Recommendation/?$', search, {"model":"Recommendation"}),
    (r'^search/SeedRecommendation/?$', search, {
		"model":"Recommendation",
		"params":{"seeded":True}
		}),

  #Users and Labs
   #Authentication
    (r'^user_logout/$', "DRP.views.user.user_logout"), #Log Out
    (r'^user_login/$', "DRP.views.user.user_login"), #Log In
    (r'^login/?$', info_page, {"page":"login_form"}),
   #Registration
    (r'^user_license_agreement/$', "DRP.views.license_agreement.get_user_license_agreement"),
    (r'^update_user_license_agreement/$', "DRP.views.license_agreement.update_user_license_agreement"),
    (r'^registration_prompt/$', "DRP.views.registration.registration_prompt"), #Redirects to correct registration choice.
    (r'^lab_registration/$', "DRP.views.registration.lab_registration"), #Create Lab ###INACTIVE
    (r'^user_registration/$', "DRP.views.registration.user_registration"), #Create User
    (r'^user_update/$', "DRP.views.user.user_update"), #Create User
    (r'^change_password/$', "DRP.views.user.change_password"), #Change Password

 #Enable the admin:
    (r'^admin/', include(admin.site.urls)),
)
