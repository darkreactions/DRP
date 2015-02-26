from django.conf.urls import patterns, url, include
from django.contrib import admin

#TODO: Abstract these out/remove these.
from core_views import upload_prompt, upload_CSV, data_form, check_compound, compound_guide_form, compound_guide_entry, edit_CG_entry, change_Recommendation, assign_user_to_rec, edit_recommendation

# Uncomment the next two lines to enable the admin: ###C
admin.autodiscover()
handler500 = 'DRP.views.errors.display_500_error'
handler404 = 'DRP.views.errors.display_404_error'
page = "(?P<page_request>\d+)"
urlpatterns = patterns('',
 #Individual Pages.
    (r'^$', "DRP.views.general.page", {"template":"home"}),
    (r'^home/?$', "DRP.views.general.page", {"template":"home"}),
    (r'^papers/?$', "DRP.views.general.page", {"template":"papers"}),
    (r'^about/?$', "DRP.views.general.page", {"template": "about"}),
    (r'^contact/?$', "DRP.views.general.page", {"template":"contact"}),
    (r'^login/?$', "DRP.views.general.page", {"template":"login_form"}),

    (r'^explore/?$', "DRP.views.general.locked_page", {"template":"explore"}),

    (r'^contact_form/?$', "DRP.views.contact.contact_form"),
  #Dashboard
    url(r'^dash(?:board)?/?$', "DRP.views.dashboard.get_dashboard"),
    url(r'^get_class_stats/(?P<classes>[24])/?$', "DRP.views.dashboard.get_class_stats_json"),

    url(r'^graphs?/?$', "DRP.views.graph.graph"),
    url(r'^graphs?/test[s_]?(?:sizes?)?/?$',"DRP.views.graph.graph",{"base":"test"}),
    url(r'^graphs?/train[s_]?(?:sizes?)?/?$',"DRP.views.graph.graph",{"base":"train"}),
    url(r'^graphs?/(?:dates?|times?)/?$', "DRP.views.graph.graph", {"base":"time"}),

  #Database
    (r'^data(?:base)?/?$', "DRP.views.database.database"), #Encompassing data view.
    (r'^data(?:base)?/'+page+'/?$', "DRP.views.database.database"), #Encompassing data view.

  # Searching
    (r'^search/$', "DRP.views.search.search", {"model":"Data"}),
    (r'^search/Data/?$', "DRP.views.search.search", {"model":"Data"}),
    (r'^search/Recommendation/?$', "DRP.views.search.search", {"model":"Recommendation"}),
    (r'^search/SeedRecommendation/?$', "DRP.views.search.search", {
      "model":"Recommendation",
      "params":{"seeded":True}
    }),


   #Upload/Download database.
    (r'^upload_prompt/$', upload_prompt),
    (r'^upload_data/$', upload_CSV),
    (r'^download_prompt/$', "DRP.views.download.download_prompt"),
    (r'^download_data/$', "DRP.views.download.download_CSV"),
   #Modify Data
    (r'^change_Data/$', "DRP.views.data_editing.change_Data"), #[JSON] Edit Data Entry
    (r'^delete_Data/$', "DRP.views.data_editing.delete_Data"), #[JSON] Delete Data Entries
    (r'^add_reactant/$', "DRP.views.data_editing.add_reactant"), #[JSON] Add Reactant Group
    (r'^delete_reactant/$', "DRP.views.data_editing.delete_reactant"), #[JSON] Delete Reactant Group
    (r'^data_form/$', data_form), #Form for adding new data.
   #Add Compound Entries
    (r'^check_compound/$', check_compound), #Form for adding new CG abbreviations.
    (r'^compound_guide/$', "DRP.views.compound_guide.compound_guide"), #Form for adding new CG abbreviations.
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

    (r'^seeds?(?:_recommend(?:ations?)?)?/?$', "DRP.views.seed_recommend.seed_recommend"),
    (r'^seeds?(?:_recommend(?:ations?)?)?/'+page+'/?$', "DRP.views.seed_recommend.seed_recommend"),

    (r'^check_seed_oven/$', "DRP.views.seed_recommend.check_seed_worker_cache"),

    (r'^rec(?:ommend(?:ation)?)?s?/?$', "DRP.views.database.database", {"model":"recommendations"}),
    (r'^rec(?:ommend(?:ation)?)?s?/'+page+'/?$',"DRP.views.database.database", {"model":"recommendations"}),

    (r'^saved/?$', "DRP.views.database.database", {"model":"saved"}),
    (r'^saved/'+page+'/?$',"DRP.views.database.database", {"model":"saved"}),

    (r'^change_Recommendation/$', change_Recommendation), #[JSON] Edit Rec Entry
    (r'^assign_user/$', assign_user_to_rec),
    (r'^show_recommendation/$', edit_recommendation, {"action":"show"}),
    (r'^hide_recommendation/$', edit_recommendation, {"action":"hide"}),
    (r'^save_recommendation/$', edit_recommendation, {"action":"save"}),
    (r'^unsave_recommendation/$', edit_recommendation, {"action":"unsave"}),
    (r'^sensical_recommendation/$', edit_recommendation, {"action":"sense"}),
    (r'^nonsensical_recommendation/$', edit_recommendation, {"action":"nonsense"}),

  #Users and Labs
   #Authentication
    (r'^user_logout/$', "DRP.views.user.user_logout"), #Log Out
    (r'^user_login/$', "DRP.views.user.user_login"), #Log In
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
