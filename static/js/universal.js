$(document).ready(function() {
//######################################################################
selectedData = Array();
currentQuery = Array();

CGEntries = undefined;
CGAbbrevs = undefined;

//############   Sorting   #############################################
//Sort from greatest to least.
window.sortNumbers = function(smallNum, bigNum) {
 return smallNum - bigNum;
}

window.sortNumbersReverse = function(smallNum, bigNum) {
 return bigNum - smallNum;
}

//############   Formatting   ##########################################
window.make_name_verbose = function(string) {
 var verbose_name = "";
 var parts = string.split("_");
 for (var i in parts){
  verbose_name += parts[i][0].toUpperCase()+parts[i].slice(1)+" ";
 }
 verbose_name = verbose_name.slice(0,-1);
 return verbose_name;
}

//Refresh the data container classes. (NOTE: Does not perform server request for new data.)
window.restyleData = function() {
 //Fade out units if the amount is also faded out.
 $(".type_notes").each(function() {
  if ($(this).is(":empty")) {
   $(this).addClass("opaqueDatum");
  }
 });

 //Keep selected data highlighted even if on page changes.
 $(".dataEntry").each(function() {
  dataID = $(this).find(".type_ref").html().trim();
  if (selectedData.indexOf(dataID) != -1) {
   $(this).parent().addClass("dataSelected");
  }
 });
}
//############   Tooltips   ############################################
//Apply custom tooltips to applicable data.
$(document).tooltip({
  position: {
  my: "center top+20",
 },
 track: true,
 content: function() {
  return $(this).attr("title")
  }
 });

//Disable the tooltips when a text-input is selected.
tooltips_disabled = false;
var exemptList = ["searchValue"];
$(document).on("focusin", "input[type=text]", function() {
 //if ((!tooltips_disabled) && exemptList.indexOf($(this).attr("id")) == -1) {
 if (exemptList.indexOf($(this).attr("id")) == -1) {
  tooltips_disabled = true
  $(document).tooltip("disable");
 }
});

$(document).on("focusout", "input", function() {
 tooltips_disabled = false
 $(document).tooltip("enable");
});

//############  Side Container:  ####################################
function toggleSideContainer() {
 if ($("#sidePanel").css("width")!="0px"){
  $("#mainPanel").css("width","100%"); 
  $("#sidePanel").css("width","0%"); 
 } else {
  $("#mainPanel").animate({"width":"50%"}, 750); 
  $("#sidePanel").animate({"width":"50%"}, 750); 
 }
}

//############   Ribbons:   #############################################
window.showRibbon = function(message, color, location, timeout) {
 //Assume that the ribbon should time out.
 timeout = timeout !== undefined ? timeout : true
 //Remove any extra
 $(".ribbonMessage").remove();

 $(location).append(
  "<div class=\"ribbonMessage\" style=\"background-color:" + color +
   ";\">" + message + "</div>"
 );
 if (timeout) {
  setTimeout(function() {
   $(location).children(".ribbonMessage").fadeOut(1000, function()
   {
    $(".ribbonMessage").remove();
   });

  },500+15*message.length);
 }
}

//############  Form Interactions:  ####################################
$(document).on("submit", ".infoForm", function() {
 var form = $(this); //Keep a reference to the form inside the POST request.
 $.post($(form).attr("action"), $(form).serialize(), function(response) {
  //Recreate the popup window with the server response.
  $("#popupContainer_inner").html(response);
  $(".subPopup").draggable();

  //Reset the autocomplete if necessary.
  $(".autocomplete_reactant").autocomplete({
   source: CGAbbrevs,
   messages: {
    noResults: "",
    results: function() {}
   }
  });

  //Show the ribbon message if applicable.
  if ($(".successActivator").length) {
   showRibbon("Data added!", "green", "#popupContainer_inner");
   $(".successActivator").remove();
   return false;
  }

  //Reload the page if applicable.
  if ($(".reloadActivator").length) {
   window.location.reload(true);
  }
  $(".loadingWheel").remove();


 });
 return false; //Do not continue or else the form will post again.
});

$(document).on("click", ".form_button[type=submit]", function() {
 $(".loadingWheel").remove();
 $(this).hide();
 $(this).parent().append("<div class=\"loadingWheel\"></div>");
});

//############ User Authentication: ####################################
$("#userLogOut").click( function() {
 //Send the log-out signal.
 $.get("/user_logout/", function() {
  //Reload the screen to verify log-out.
  window.location.reload(true);
  });
});

//################   Search   ##########################################
$(document).on("click", ".PT_element", function() {
 $(this).toggleClass("PT_selected");

 var element = $(this).html().trim();
 var indexOfElement = PT_selected.indexOf(element);

 if (indexOfElement == -1){
  PT_selected.push(element);
 } else {
  PT_selected.splice(indexOfElement,1);
 }
});

function sendSearchQuery(currentQuery) {
 JSONQuery = JSON.stringify({"currentQuery":currentQuery,"page":"1"});

 $.ajax({
  url:"/search/",
  method:"post",
  data: {"body":JSONQuery},
  traditional: true,
  success: function(response) {
   $("#mainPanel").html(response)
   restyleData();
  }
 });
}

//Back button tooltip
$(document).on("mouseover", ".search_backButton", function() {
 if (currentQuery.length){
  //Get the previously used filters.
  var filter_string = "Filters:"
  for (var i in currentQuery.slice(0,-1)) {
   filter_string += "<br/>"+(parseInt(i)+1)+".) "+make_name_verbose(currentQuery[i]["field"])+": " + currentQuery[i]["value"]
  }

  //Add the new filter.
  filter_string += "<br/><div class=\"search_backText\">"+parseInt(currentQuery.length)+".) "+make_name_verbose(currentQuery[currentQuery.length-1]["field"])+": "+currentQuery[currentQuery.length-1]["value"]+"</div>"

  $(this).attr("title", filter_string);
 } else {
  $(this).attr("title", "Remove the last filter.");
 }
});

//Filter button tooltip
function get_atom_query() {
 var array_atom_query = Array();

 $(".PT_selected").each(function() {
  array_atom_query.push($(this).html());
 });

 var current_atom_query = array_atom_query.join(", ");

 if (array_atom_query.length>2) {
  var i = current_atom_query.lastIndexOf(" ");
  current_atom_query = [current_atom_query.slice(0, i), " " + $(".radioSearch:checked").val(), current_atom_query.slice(i)].join(" ");
 } else if (array_atom_query.length==2) {
  current_atom_query = current_atom_query.replace(",", " " + $(".radioSearch:checked").val());
 }

 return current_atom_query.replace(/ +/g, " ")
}

$(document).on("mouseover", ".search_filterButton", function() {
 //Get the previously used filters.
 var filter_string = "Filters:"
 for (var i in currentQuery) {
  if (currentQuery[i]["field"]=="atoms"){
   filter_string += "<br/>"+(parseInt(i)+1)+".) "+make_name_verbose(currentQuery[i]["field"])+": "+currentQuery[i]["value"];
  } else {
   filter_string += "<br/>"+(parseInt(i)+1)+".) "+make_name_verbose(currentQuery[i]["field"])+": \""+currentQuery[i]["value"]+"\" ("+currentQuery[i]["match"]+")</div>";
  }
 }
 
 //If "Atoms" search is active.
 if ($(".ui-state-active").children().html()=="Atoms" && $(".PT_selected").length>0) {
  current_atom_query = get_atom_query();

  //Add the new filter.
  filter_string += "<br/><div class=\"search_filterText\">"+(parseInt(currentQuery.length)+1)+".) "+"Atoms: "+current_atom_query+"</div>"

 } else if ($(".ui-state-active").children().html()=="Fields" && $("#searchValue").val()){
   field = $(".dropDownMenu[name=field] option:selected").val();
   match = $(".dropDownMenu[name=match] option:selected").val();
   value = $("#searchValue").val();
   //Add the new filter.
   filter_string += "<br/><div class=\"search_filterText\">"+(parseInt(currentQuery.length)+1)+".) "+make_name_verbose(field)+": \""+value+"\" ("+match+")</div>";
 } else {
  filter_string = "Enter a query!"; 
 }
 //Apply the title.
 $(this).attr("title", filter_string);
});

$(document).on("click", ".search_filterButton", function() {
 if ($(".dataGroup").length != 0) {
  //If "Atoms" search is active.
  if ($(".ui-state-active").children().html()=="Atoms" && $(".PT_selected").length>0) {
   field = "atoms";
   value = get_atom_query().replace(/,/g,"");
   //Make sure the query was not already searched.
   if (currentQuery){
    for (var i in currentQuery) {
     if (currentQuery[i]["field"] == field && currentQuery[i]["value"] == value) {
      showRibbon("Already queried!", "#FFC87C","#sidePanel", true);
      return false //Don't continue if query is already present.
     }
    }
   }
   showRibbon("Searching!", "#99FF5E","#sidePanel", true);

   currentQuery.push({
    "field":field,
    "match":"exact",
    "value":value,
   });
   sendSearchQuery(currentQuery);

  //If "Fields" search is active.
  } else if ($(".ui-state-active").children().html()=="Fields" && $("#searchValue").val()){
   field = $(".dropDownMenu[name=field] option:selected").val()
   match = $(".dropDownMenu[name=match] option:selected").val()
   value = $("#searchValue").val()
   //Make sure the query was not already searched.
   if (currentQuery){
    for (var i in currentQuery) {
     if (currentQuery[i]["field"] == field && currentQuery[i]["value"] == value) {
      showRibbon("Already queried!", "#FFC87C","#sidePanel", true);
      return false //Don't continue if query is already present.
     }
    }
   }
   showRibbon("Searching!", "#99FF5E","#sidePanel", true);

   currentQuery.push({
    "field":field,
    "match":match,
    "value":value,
   });
   sendSearchQuery(currentQuery);
  } else {
   showRibbon("Nothing entered!", "#FF6870", $("#sidePanel"), true);
  }
 } else {
  showRibbon("No data to filter!", "#FF6870", $("#sidePanel"), true);
 }
});

$(document).on("click", ".search_backButton", function() {
 if (currentQuery.length){
  currentQuery.pop();
  showRibbon("Removing last filter!", "#99FF5E","#sidePanel", true);
  sendSearchQuery(currentQuery);
 } else {
  showRibbon("No filters present!", "#FF6870", $("#sidePanel"), true);
 }
});

$(document).on("click", ".search_clearButton", function() {
 currentQuery = Array();
 sendSearchQuery(currentQuery);
 showRibbon("Filters emptied", "#99FF5E","#sidePanel", true);
});

 //################   Autocompleting Search   ####################### //###
window.setAutoComplete = function(location, source) {
 $(location).autocomplete({
  source: source,
  messages: {
   noResults: "",
   results: function() {}
  }
 });
}

$(document).on("focusin", "#searchValue", function() {
 //Set autocomplete to vary with selection.
});

 //////source = source !== undefined ? source : CGAbbrevs;
 //////if (source==CGAbbrevs) {
  ////////Get the auto-complete options form the CG guide if no other source is found.
  //////if (CGAbbrevs == undefined) {
   //////$.get("/send_CG_names/", function(response) {
    //////CGEntries = response;
    //////CGAbbrevs = Array();
    //////for (var key in CGEntries) {
     //////CGAbbrevs.push(key);
    //////}
   //////})
  //////} else {
   //////$(location).autocomplete({
    //////source: CGAbbrevs,
    //////messages: {
     //////noResults: "",
     //////results: function() {}
    //////}
   //////});
  //////}
 //////} else {
  //////$(location).autocomplete({
   //////source: CGAbbrevs,
   //////messages: {
    //////noResults: "",
    //////results: function() {}
   //////}
  //////});
 //////}
//////}//###


//############ Popup Management: #######################################
window.createPopupConfirmation = function(message) {
 $("body").append("<div class=popupConfirmation title=Confirm:>"+
  message+"</div>");
}

//Activate specified popup:
$(document).on("click", ".popupActivator", function(event) {
 //Display a loading message if the request takes a visible amount of time.
 $("#popupContainer_inner").html("<center>Loading. Please wait.</center>");

 //Load the activator CSS.
 activatorID = $(this).attr("id");

 //"Forward" specific triggers to other popups.
 var specificRef = "";
 if (activatorID == "leftMenu_addNew_copy") {
  activatorID = "leftMenu_addNew"
  specificDatum =  $(this).siblings(".dataEntry").children(".type_ref")//Mark spaces in references by "+"
  specificRef = $(specificDatum).html().trim().replace(" ","+");
  event.stopPropagation();
 }

 $("#popupContainer").attr("for", activatorID);
 switch (activatorID) {
  case "leftMenu_addNew":
   //Send the request to the server
   $.get("/data_form/"+specificRef, function(response) {
    $("#popupContainer_inner").html(response);

    //Get the auto-complete options form the CG guide. ###SOME ERROR HERE, CASEY
    if (CGAbbrevs == undefined) {
     $.get("/send_CG_names/", function(response) {
      CGEntries = response;
      CGAbbrevs = Array();
      for (var key in CGEntries) {
       CGAbbrevs.push(key);
      }
      $(".autocomplete_reactant").autocomplete({
       source: CGAbbrevs,
       messages: {
        noResults: "",
        results: function() {}
       }
      });
     })
    } else {
     $(".autocomplete_reactant").autocomplete({
      source: CGAbbrevs,
      messages: {
       noResults: "",
       results: function() {}
      }
     });
    }
   });
   break;
  case "leftMenu_uploadCSV":
   $.get("/upload_CSV/", function(response) {
    $("#popupContainer_inner").html(response);
   });
   break;
  case "leftMenu_downloadCSV":
   $.get("/download_CSV/", function(response) {
    $("#popupContainer_inner").html(response);
   });
   break;
  case "searchButton":
   PT_selected = Array();
   $.get("/search/", function(response) {
    $("#sidePanel_inner").html(response);
    toggleSideContainer();
    $("#tabs").tabs({active: 1});
   });
   return false; //TODO: Separate this into a "sidePanel" activator
   break;
  case "userLogin":
   $.get("/user_login/", function(response) {
    $("#popupContainer_inner").html(response);
   });
   break;
  case "changePassword":
   $.get("/change_password/", function(response) {
    $("#popupContainer_inner").html(response);
   });
   break;
  case "registrationPrompt": //###REPLACE WITH PRELOADED DROP-DOWN MENU?
   $.get("/registration_prompt/", function(response) {
    $("#popupContainer_inner").html(response);
   });
   break;
  case "userUpdate":
   $.get("/user_update/", function(response) {
    $("#popupContainer_inner").html(response);
   });
   break;
  case "userRegistration":
   $.get("/user_registration/", function(response) {
    $("#popupContainer_inner").html(response);
   });
   break;
  case "labRegistration":
   $.get("/lab_registration/", function(response) {
    $("#popupContainer_inner").html(response);
   });
   break;
  case "CG_activatorButton":
   $.get("/compound_guide_form/", function(response) {
    $("#popupContainer_inner").html(response);
    $(".CG_saveButton").remove()
    //Erase the current JSON object and selected CG entries.
    CGEntries = undefined;
    CGAbbrevs = undefined;
    CGSelected = Array();
   });
   break;
  default:
   return false; //If a popup is not recognized, don't load anything.
 }
 $("#popupGlobal").fadeIn(300);
});

// Fade the popup when the mask is clicked.
$(document).on("click", "#mask", function() {
 if ($(".ribbonMessage").length==0 && $(".maskBlockFade").length==0){
  $("#popupGlobal").fadeOut("fast");

  //Remove any extra additions the popup may have populated.
  $(".CG_saveButton").remove()

  //Reload the screen if requested.
  if ($(".reloadActivator").length) {
   window.location.reload(true);
   $(".reloadActivator").remove();
  };
 }
});

function darkenMask() {
 $("#mask").addClass("darkenedMask");
}

function revertMask() {
 $("#mask").removeClass("darkenedMask");
 $(".darkenedMask").remove();
}


// Cancel any masks or popups that are covering the screen. 
$(document).on("click", ".clearScreenButton", function() {
  $("#popupGlobal").fadeOut("fast");

  //Remove any extra additions the popup may have populated.
  $(".CG_saveButton").remove()
  $(".closeButton").show();
  $("#mask").removeClass("darkenedMask");
  $(".maskBlockFade").remove();
});

//Close popups on close-button click.
$(document).on("click", ".closeButton", function() {
 if ($(".ribbonMessage").length==0 && $(".maskBlockFade").length==0){
  //Close the global container if the main container is closed.
  CGSelected = Array();
  if ($(this).hasClass("refreshOnDie")){
   window.location.reload(true);
  } else if ($(this).parent().attr("id")=="popupContainer") {
   $(this).parent().parent().fadeOut("fast");
  } else {
   $(this).parent().fadeOut("fast");
  }
 }
});

//Shrink the slide container on contractButton clicks.
$(document).on("click", ".contractButton", function() {
 $("#sidePanelWarning").hide();
 toggleSideContainer();
});

//Detect any activators that might be loaded.
if ($(".loadActivator").length){
 var activatorID = $(".loadActivator").attr("id");
 $("#popupContainer").attr("for", activatorID);
 switch (activatorID) {
  case "userLicenseAgreement":
   darkenMask();
   $.get("/user_license_agreement/", function(response) {
    $("#popupContainer_inner").html(response);
    $(".closeButton").hide(); //Remove unnecessary close buttons.
   });
   break;
 }; 
 $("#popupGlobal").fadeIn(300);
}

//######################################################################
});








