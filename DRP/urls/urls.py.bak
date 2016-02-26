from django.conf.urls import url, include
from django.contrib import admin

#TODO: Abstract these out/remove these.
from core_views import upload_prompt, upload_CSV, data_form, check_compound, compound_guide_form, compound_guide_entry, edit_CG_entry, change_Recommendation, assign_user_to_rec, edit_recommendation

# Uncomment the next two lines to enable the admin: ###C
admin.autodiscover()
handler500 = 'DRP.views.errors.display_500_error'
handler404 = 'DRP.views.errors.display_404_error'
page = "(?P<page_request>\d+)"
urlpatterns = [
 #Individual Pages.
    url(r'^$', DRP.views.general.page, {"template":"home"}),
    url(r'^home/?$', DRP.views.general.page, {"template":"home"}),
    url(r'^papers/?$', DRP.views.general.page, {"template":"papers"}),
    url(r'^about/?$', DRP.views.general.page, {"template": "about"}),
    url(r'^contact/?$', DRP.views.general.page, {"template":"contact"}),
    url(r'^login/?$', DRP.views.general.page, {"template":"login_form"}),
    url(r'^explore/?$', DRP.views.general.page, {"template":"explore"}),

    url(r'^contact_form/?$', DRP.views.contact.contact_form),
  #Dashboard
    url(r'^dash(?:board)?/?$', DRP.views.dashboard.get_dashboard),
    url(r'^get_class_stats/(?P<category>[24]-(?:test|train))/?$', DRP.views.dashboard.get_class_stats_json),

    url(r'^graphs?/?$', DRP.views.graph.graph),
    url(r'^graphs?/test[s_]?(?:sizes?)?/?$',DRP.views.graph.graph,{"base":"test"}),
    url(r'^graphs?/train[s_]?(?:sizes?)?/?$',DRP.views.graph.graph,{"base":"train"}),
    url(r'^graphs?/(?:dates?|times?)/?$', DRP.views.graph.graph, {"base":"time"}),

    url(r'^models?/(?P<model_id>\d+)/?$', DRP.views.stats.model_stats),

  #Database
    url(r'^data(?:base)?/?$', DRP.views.database.database), #Encompassing data view.
    url(r'^data(?:base)?/'+page+'/?$', DRP.views.database.database), #Encompassing data view.

  # Searching
    url(r'^search/$', DRP.views.search.search, {"model":"Data"}),
    url(r'^search/Data/?$', DRP.views.search.search, {"model":"Data"}),
    url(r'^search/Explore/?$', DRP.views.explore_search.search, {"model":"Data"}),
    url(r'^search/Recommendation/?$', DRP.views.search.search, {"model":"Recommendation"}),
    url(r'^search/SeedRecommendation/?$', DRP.views.search.search, {
      "model":"Recommendation",
      "params":{"seeded":True}
    }),


   #Upload/Download database.
    url(r'^upload_prompt/$', upload_prompt),
    url(r'^upload_data/$', upload_CSV),
    url(r'^download_prompt/$', DRP.views.download.download_prompt),
    url(r'^download_data/$', DRP.views.download.download_CSV),
   #Modify Data
    url(r'^change_Data/$', DRP.views.data_editing.change_Data), #[JSON] Edit Data Entry
    url(r'^delete_Data/$', DRP.views.data_editing.delete_Data), #[JSON] Delete Data Entries
    url(r'^add_reactant/$', DRP.views.data_editing.add_reactant), #[JSON] Add Reactant Group
    url(r'^delete_reactant/$', DRP.views.data_editing.delete_reactant), #[JSON] Delete Reactant Group
    url(r'^data_form/$', data_form), #Form for adding new data.
   #Add Compound Entries
    url(r'^check_compound/$', check_compound), #Form for adding new CG abbreviations.
    url(r'^compound_guide/$', DRP.views.compound_guide.compound_guide), #Form for adding new CG abbreviations.
    url(r'^compound_guide_form/$', compound_guide_form), #Form for adding new CG abbreviations.
    url(r'^compound_guide_entry/$', compound_guide_entry), #Return a single CG table entry.
    url(r'^edit_CG_entry/$', edit_CG_entry), #Edit a CG entry.
   #Validation
    url(r'^send_CG_names/$', DRP.views.jsonViews.send_CG_names), #Send the CG name_pairs for client-side validation.
  #Visualization
    url(r'^get_graph/$', DRP.views.explore_vis.get_graph_data),
    url(r'^setup_graph/$', DRP.views.explore_vis.store_graph),
  #Recommendations
    url(r'^make_seed_recommendations/$', DRP.views.seed_recommend.make_seed_recommendations),

    url(r'^seeds?(?:_recommend(?:ations?)?)?/?$', DRP.views.seed_recommend.seed_recommend),
    url(r'^seeds?(?:_recommend(?:ations?)?)?/'+page+'/?$', DRP.views.seed_recommend.seed_recommend),

    url(r'^check_seed_oven/$', DRP.views.seed_recommend.check_seed_worker_cache),

    url(r'^rec(?:ommend(?:ation)?)?s?/?$', DRP.views.database.database, {"model":"recommendations"}),
    url(r'^rec(?:ommend(?:ation)?)?s?/'+page+'/?$',DRP.views.database.database, {"model":"recommendations"}),

    url(r'^saved/?$', DRP.views.database.database, {"model":"saved"}),
    url(r'^saved/'+page+'/?$',DRP.views.database.database, {"model":"saved"}),

    url(r'^change_Recommendation/$', change_Recommendation), #[JSON] Edit Rec Entry
    url(r'^assign_user/$', assign_user_to_rec),
    url(r'^show_recommendation/$', edit_recommendation, {"action":"show"}),
    url(r'^hide_recommendation/$', edit_recommendation, {"action":"hide"}),
    url(r'^save_recommendation/$', edit_recommendation, {"action":"save"}),
    url(r'^unsave_recommendation/$', edit_recommendation, {"action":"unsave"}),
    url(r'^sensical_recommendation/$', edit_recommendation, {"action":"sense"}),
    url(r'^nonsensical_recommendation/$', edit_recommendation, {"action":"nonsense"}),

  #Users and Labs
   #Authentication
    url(r'^user_logout/$', DRP.views.user.user_logout), #Log Out
    url(r'^user_login/$', DRP.views.user.user_login), #Log In
   #Registration
    url(r'^user_license_agreement/$', DRP.views.license_agreement.get_user_license_agreement),
    url(r'^update_user_license_agreement/$', DRP.views.license_agreement.update_user_license_agreement),
    url(r'^registration_prompt/$', DRP.views.registration.registration_prompt), #Redirects to correct registration choice.
    url(r'^lab_registration/$', DRP.views.registration.lab_registration), #Create Lab ###INACTIVE
    url(r'^user_registration/$', DRP.views.registration.user_registration), #Create User
    url(r'^user_update/$', DRP.views.user.user_update), #Create User
    url(r'^change_password/$', DRP.views.user.change_password), #Change Password

 #Enable the admin:
    url(r'^admin/', include(admin.site.urls)),
)
