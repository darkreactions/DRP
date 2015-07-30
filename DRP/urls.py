from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.views.generic.base import TemplateView
from django.contrib.auth.views import login, logout
import DRP.views

#TODO: Abstract these out/remove these.
#from core_views import upload_prompt, upload_CSV, data_form, check_compound, compound_guide_form, compound_guide_entry, edit_CG_entry, change_Recommendation, assign_user_to_rec, edit_recommendation

# Uncomment the next two lines to enable the admin: ###C
admin.autodiscover()
handler500 = 'DRP.views.errors.display_500_error' #TODO: check what these views do
handler404 = 'DRP.views.errors.display_404_error' #TODO: is this the right place for this (should be set in settings.py)?
page = "(?P<page_request>\d+)"
urlpatterns = patterns('',
 #Public Pages.
    url(r'^$', TemplateView.as_view(template_name="home.html"), name='home'),
    (r'^about.html$', TemplateView.as_view(template_name="about.html")),
    url(r'^contact.html$', DRP.views.contact, name='contact'),
 #Authentication pages
    url(r'^login.html$', login, {'template_name':'login.html'}, name='login'),
    (r'^logout.html$', logout, {'next_page':'home'}),
    url(r'^register.html$', DRP.views.register, name='register'),
    url(r'^confirm.html$', DRP.views.confirm, name='confirm'),
#    (r'^license.html$', DRP.views.license)
#    (r'^join_lab.html$', DRP.views.joinLab)

#    (r'^explore.html$', "DRP.views.general.page", {"template":"explore"}),
#
#  #Dashboard
#    (r'^dashboard.html$', "DRP.views.dashboard.get_dashboard"),
#  #Database
#    (r'^database/$', "DRP.views.database.database"), #Encompassing data view.
#    (r'^database/'+page+'/?$', "DRP.views.database.database"), #Encompassing data view.
#
#  # Searching
#    (r'^search/$', "DRP.views.search.search", {"model":"Data"}),
#    (r'^search/data/?$', "DRP.views.search.search", {"model":"Data"}),
#    (r'^search/explore/?$', "DRP.views.explore_search.search", {"model":"Data"}),
#    (r'^search/recommendation/?$', "DRP.views.search.search", {"model":"Recommendation"}),
#    (r'^search/seed_recommendation/?$', "DRP.views.search.search", {
#      "model":"Recommendation",
#      "params":{"seeded":True}
#    }),
#
#
#   #Upload/Download database.
#    (r'^upload_data/$', upload_CSV),
#    (r'^download_prompt/$', "DRP.views.download.download_prompt"),
#    (r'^download_data/$', "DRP.views.download.download_CSV"),
   #Modify Data
#    (r'^change_Data/$', "DRP.views.data_editing.change_Data"), #[JSON] Edit Data Entry
#    (r'^delete_Data/$', "DRP.views.data_editing.delete_Data"), #[JSON] Delete Data Entries
#    (r'^add_reactant/$', "DRP.views.data_editing.add_reactant"), #[JSON] Add Reactant Group
#    (r'^delete_reactant/$', "DRP.views.data_editing.delete_reactant"), #[JSON] Delete Reactant Group
#    (r'^data_form/$', data_form), #Form for adding new data.
   #Add Compound Entries
#    (r'^check_compound/$', check_compound), #Form for adding new CG abbreviations.
#    (r'^compound_guide/$', "DRP.views.compound_guide.compound_guide"), #Form for adding new CG abbreviations.
#    (r'^compound_guide_form/$', compound_guide_form), #Form for adding new CG abbreviations.
#    (r'^compound_guide_entry/$', compound_guide_entry), #Return a single CG table entry.
#    (r'^edit_CG_entry/$', edit_CG_entry), #Edit a CG entry.
   #Validation
#    (r'^send_CG_names/$', "DRP.views.jsonViews.send_CG_names"), #Send the CG name_pairs for client-side validation.
  #Visualization
#    (r'^get_graph/$', "DRP.views.explore_vis.get_graph_data"),
#    (r'^setup_graph/$', "DRP.views.explore_vis.store_graph"),
  #Recommendations
#    (r'^make_seed_recommendations/$', "DRP.views.seed_recommend.make_seed_recommendations"),

 #   (r'^seeds?(?:_recommend(?:ations?)?)?/?$', "DRP.views.seed_recommend.seed_recommend"),
#    (r'^seeds?(?:_recommend(?:ations?)?)?/'+page+'/?$', "DRP.views.seed_recommend.seed_recommend"),

 #   (r'^check_seed_oven/$', "DRP.views.seed_recommend.check_seed_worker_cache"),

 #   (r'^rec(?:ommend(?:ation)?)?s?/?$', "DRP.views.database.database", {"model":"recommendations"}),
 #   (r'^rec(?:ommend(?:ation)?)?s?/'+page+'/?$',"DRP.views.database.database", {"model":"recommendations"}),

  #  (r'^saved/?$', "DRP.views.database.database", {"model":"saved"}),
#    (r'^saved/'+page+'/?$',"DRP.views.database.database", {"model":"saved"}),

#    (r'^change_Recommendation/$', change_Recommendation), #[JSON] Edit Rec Entry
#    (r'^assign_user/$', assign_user_to_rec),
#    (r'^show_recommendation/$', edit_recommendation, {"action":"show"}),
#    (r'^hide_recommendation/$', edit_recommendation, {"action":"hide"}),
#    (r'^save_recommendation/$', edit_recommendation, {"action":"save"}),
#    (r'^unsave_recommendation/$', edit_recommendation, {"action":"unsave"}),
#    (r'^sensical_recommendation/$', edit_recommendation, {"action":"sense"}),
#    (r'^nonsensical_recommendation/$', edit_recommendation, {"action":"nonsense"}),

  #Users and Labs
   #Authentication
#    (r'^user_logout/$', "DRP.views.user.user_logout"), #Log Out
#    (r'^user_login/$', "DRP.views.user.user_login"), #Log In
   #Registration
#    (r'^user_license_agreement/$', "DRP.views.license_agreement.get_user_license_agreement"),
#    (r'^update_user_license_agreement/$', "DRP.views.license_agreement.update_user_license_agreement"),
#    (r'^registration_prompt/$', "DRP.views.registration.registration_prompt"), #Redirects to correct registration choice.
#    (r'^lab_registration/$', "DRP.views.registration.lab_registration"), #Create Lab ###INACTIVE
#    (r'^user_registration/$', "DRP.views.registration.user_registration"), #Create User
#    (r'^user_update/$', "DRP.views.user.user_update"), #Create User
#    (r'^change_password/$', "DRP.views.user.change_password"), #Change Password

 #Enable the admin:
    (r'^admin/', include(admin.site.urls)),
)
