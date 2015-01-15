//Pre-load Variable Setup:
goodColor = "#99FF5E";
neutralColor = "#FFC87C";
badColor = "#FF6870";

$(document).ready(function() {
//######################################################################
selectedData = Array();
currentQuery = Array();

CGEntries = undefined;
CGAbbrevs = undefined;

//############   Dependent Functions  ##################################
function adaptSize(element) {
 var numChars = $(element).val().length;
 var maxWidth = $(element).css("max-width") !== "none" ? parseInt($(element).css("max-width")) : 400 ;
 var proposedWidth = ((numChars*6)+35)
 if (proposedWidth < maxWidth) {
  var newWidth = proposedWidth;
 } else {
  var newWidth = maxWidth;
 }
 $(element).animate({
  "width": newWidth,
 }, 50);
};

function formatText(text) {
 return text.replace(/_{(.*?)}/g, function(r,m){return "<span class=\"subscript\">"+m+"</span>"});
}
window.applyTickingDots = function(location) {
  var interval = setInterval(function(){
    content = $(location).html();
    dotCount = content.split(".");
    if (dotCount.length >= 4){
      $(location).html(".");
    } else {
      $(location).append(" .");
    }
  }, 1*1000);

  //Bind the loading animation to this specific entry.
  $(location).data("timer", interval);

  return interval; //Return the interval so that it may be cleared.
}

function clearLoadingWheels(container){
  //Stop the timers in the soon-to-be-removed .loadingWheels.
  clearInterval($(container).find(".loadingWheel").data("timer"));

  //Remove any pre-existing .loadingWheels in the container.
  $(container).find(".loadingWheel").remove();
}

window.applyLoadingWheel = function(location){
  //Remove old .loadingWheels if they exist.
  clearLoadingWheels(location);

  //Create a new .loadingWheel in the container and return the load timer.
  var newLoad = $("<div class=\"loadingWheel\">.</div>").appendTo(location);
  return applyTickingDots(newLoad);
}

//Get edit-by-menu options from editChoices.json.
function getOptions(field) {
 switch (field) {
  case ("unit"):
   return  editChoices["unitChoices"]
  case ("compound_type"):
   return  editChoices["typeChoices"]
  case ("recommended"):
   return  editChoices["boolChoices"]
  case ("outcome"):
   return  editChoices["outcomeChoices"]
  case ("purity"):
   return  editChoices["purityChoices"]
  case ("slow_cool"):case("leak"): //ie, slow_cool and leak are the same case.
   return  editChoices["boolChoices"]
  default:
   return ["No options found."]
 }
}

function reloadIfNeeded(){
 //Reload the page immediately if applicable.
 if ($(".reloadImmediately").length) {
  if ($(".reloadImmediately").attr("name")==""){
    window.location.reload(true);
  } else {
    window.location.href = $(".reloadImmediately").attr("redirect");
  }
 }
}

window.refreshOnMaskFade = function() {
 $("body").append("<div class=\"reloadActivator\"></div>")
}

function getFormErrors(form) {
 var formID = $(form).attr("action");
 var errors = [];
 switch(formID) {
  case ("/contact_form/"):
   //Verify the email.
   var emailRegExp = new RegExp(/^\S+@\S+\.\S+$/);
   var email = $(form).find("#email").val();
   if (! emailRegExp.test(email)){
     errors.push("Email invalid!");
   }

   //Verify that content and subject exist.
   var content = $(form).find("#content").val();
   if (content==""){
     errors.push("No message written!");
   }

   var subject = $(form).find("#subject").val();
   if (subject==""){
     errors.push("No subject input!");
   }

  default: //If not listed, then assume no validation is needed.
   return errors;
 }
}

function loadPageFilters(){
  $(".filterLoader").each(function(i){
    currentQuery.push({
      "field":$(this).attr("field"),
      "match":"exact",
      "value":$(this).attr("value"),
    });
  });
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

 try {
  //Keep selected data highlighted even if on page changes.
  $(".dataGroup").each(function() {
   dataID = $(this).attr("pid");
   if (selectedData.indexOf(dataID) != -1) {
    $(this).addClass("dataSelected");
   }
  });
 } catch(err){}
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
  },
 fade: 400
 });

//Disable the tooltips when a text-input is selected.
var exemptList = ["searchValue"]
tooltips_disabled = false;
$(document).on("focusin", "input[type=text]", function() {
 if ((!tooltips_disabled) && exemptList.indexOf($(this).attr("id")) == -1) {
  tooltips_disabled = true
  $(document).tooltip("disable");
 }
});

$(document).on("focusout", "input", function() {
 tooltips_disabled = false
 $(document).tooltip("enable");
});

//Reactant Tooltips.
$(document).on("mouseover", ".type_reactant", function() {
 var specificDiv = $(this);
 //If the reactant is empty, don't complain about the compound guide.
 if ($(specificDiv).is(":empty") || $(".editField").length){
  $(specificDiv).attr("title", "")
 } else if (CGEntries == undefined) {
  $(specificDiv).attr("title", "Loading!")
  $.get("/send_CG_names/", function(response) {
   CGEntries = response;
   var compound = CGEntries[$(specificDiv).html().trim()] || "Compound not in guide!"
   $(".ui-tooltip").html(compound);
   });
 } else {
  var compound = CGEntries[$(this).html().trim()] || "Compound not in guide!"
  $(this).attr("title", compound)
 }
});
//############  Side Container:  ####################################
function toggleSideContainer() {
 $("#sidePanel").toggle();
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

//############  Mask Interactions:  ####################################
function darkenMask() {
 $("#mask").addClass("darkenedMask");
 $("#mask").addClass("maskBlockFade");
 $(".closeButton").hide();
}

function revertMask() {
 $("#mask").removeClass("darkenedMask");
 $("#mask").removeClass("maskBlockFade");
 $(".closeButton").show();
}

// Fade the popup when the mask is clicked.
$(document).on("click", "#mask", function() {
 if ($(".ribbonMessage").length==0 && $(".maskBlockFade").length==0){
  $("#popupGlobal").fadeOut("fast");

  //Remove any extra additions the popup may have populated.
  $(".CG_saveButton").remove()

  //Remove any edits that are pending.
  $(".editField").each(function() {
   var oldVal = $(this).attr("oldVal");
   $(this).parent().html(oldVal);
   });

  //Reload the screen if requested.
  if ($(".reloadActivator").length) {
   window.location.reload(true);
   $(".reloadActivator").remove();
  };
 }
});

//############  Button Propagation:  ####################################
window.addDataSpecificButton = function(data, buttonID, image, title, classes) {
 classes = classes !== undefined ? classes : "";

 var buttonDiv = "<div id=\""+buttonID+"\""
 buttonDiv += "class=\"dataSpecificButton "+classes+"\"";
 buttonDiv += "style=\"background-image: url(";
 buttonDiv += STATIC_URL+"/icons/"+image+");\""
 buttonDiv += "title=\""+title +"\"></div>";
 $(data).find(".dataSpecificButtonContainer").append(buttonDiv);
}

$(document).on("mouseleave", ".dataGroup", function() {
 //Eliminate any duplication buttons on when the mouse isn't on a group.
 $(".dataSpecificButton").remove();
});

$(document).on("mouseover", ".dataGroup", function() {
 if ($(".dataSpecificButton").length == 0 ){

  //DATA-SPECIFIC BUTTONS # # # # # # # # # # # # # # # #
  //Add the copy button.
  if ($(this).attr("class").indexOf("copyable") >= 0) {

   addDataSpecificButton(this, "leftMenu_addNew", "add.png",
    "Copy this reaction to the data form.",
    "popupActivator duplicateSpecificDataButton");
   addDataSpecificButton(this, "makeSeedRecommendations", "seed.gif",
    "Generate recommendations based on this datum.",
    "");

  }

  //RECOMMENDATION-SPECIFIC BUTTONS # # # # # # # # # # # # # # # #
  if ($(this).attr("class").indexOf("recommendation") >= 0) {
   //Add the (un)save buttons.
   if ($(this).attr("class").indexOf("savedRecommendation") < 0) {
    addDataSpecificButton(this, "saveRecommendation", "save.png",
     "Save this recommendation.")
   } else {
    addDataSpecificButton(this, "unsaveRecommendation", "delete.png",
     "Unsave this recommendation.")
   }

   //Add the (non)sensical buttons.
   if ($(this).attr("class").indexOf("badRecommendation") < 0) {
    addDataSpecificButton(this, "nonsensicalRecommendation", "nonsense.png",
     "Mark this recommendation as nonsensical.")
   } else {
    addDataSpecificButton(this, "sensicalRecommendation", "check.png",
     "Mark this recommendation as sensical.")
   }

   //Add the show/hide button.
   if ($(this).attr("class").indexOf("hiddenRecommendation") < 0) {
    addDataSpecificButton(this, "hideRecommendation", "hide.png",
     "Hide this recommendation.");
   } else {
    addDataSpecificButton(this, "showRecommendation", "show.png",
     "Show this recommendation again.");
   }

   //Add the transfer-to-database buttons.
   if ($(this).attr("class").indexOf("transferable") >= 0) {
    addDataSpecificButton(this, "transferRecommendation", "add.png",
     "Complete this entry.",
     "popupActivator");
   }

  }

 }
});

$(document).on("click", ".dataSpecificButton", function(event) {
 var dataGroup = $(this).closest(".dataGroup");
 var buttonID = $(this).attr("id");

 try {
  var pid = $(dataGroup).attr("pid");
 } catch(err) {
  return false;
 }

 //Get the action of the button.
 switch (buttonID){
  case ("makeSeedRecommendations"):
   var url="/make_seed_recommendations/";
   var JSON = {"pid":pid}
   break;
  case ("saveRecommendation"):
   var url="/save_recommendation/";
   var JSON = {"pid":pid}
   break;
  case ("unsaveRecommendation"):
   var url="/unsave_recommendation/";
   var JSON = {"pid":pid}
   break;
  case ("showRecommendation"):
   var url="/show_recommendation/";
   var JSON = {"pid":pid}
   break;
  case ("hideRecommendation"):
   var url="/hide_recommendation/";
   var JSON = {"pid":pid}
   break;
  case ("sensicalRecommendation"):
   var url="/sensical_recommendation/";
   var JSON = {"pid":pid}
   break;
  case ("nonsensicalRecommendation"):
   var url="/nonsensical_recommendation/";
   var JSON = {"pid":pid}
   break;
  default:
   return false;
 }


 //Send the request and do something with the response
 $.post(url, JSON, function(response) {
  if (response=="0"){
   switch (buttonID) {
    case("makeSeedRecommendations"):
     var comment = "Making recommendations based on seed!";
     break;
    case("saveRecommendation"):
     $(dataGroup).addClass("savedRecommendation");
     var comment = "Saved!";
     break;
    case("unsaveRecommendation"):
     $(dataGroup).removeClass("savedRecommendation");
     var comment = "Unsaved!";
     break;
    case("hideRecommendation"):
     $(dataGroup).addClass("hiddenRecommendation");
     var comment = "Hidden!";
     break;
    case("showRecommendation"):
     $(dataGroup).removeClass("hiddenRecommendation");
     var comment = "Shown!";
     break;
    case("sensicalRecommendation"):
     $(dataGroup).removeClass("badRecommendation");
     var comment = "Marked as sensical!";
     break;
    case("nonsensicalRecommendation"):
     $(dataGroup).addClass("badRecommendation");
     var comment = "Marked as nonsense!";
     break;
   }

   //Refresh the mouse buttons.
   $(dataGroup).children(".dataSpecificButtonContainer").empty();
   $(dataGroup).trigger("mouseover");
   showRibbon(comment, goodColor, "#mainPanel");
  } else {
   var failureMessage;
   switch (buttonID) {

    case("makeSeedRecommendations"):
     if (response=="2"){
      failureMessage = "Still working on the last batch of recommendations!";
     } else {
      failureMessage = "Could not make recommendations from seed!";
     }
     break;

    default:
     failureMessage = "Edit failed!";
   }
   showRibbon(failureMessage, badColor, "#mainPanel");
  }

 });

 event.stopPropagation();

});

//Cancel Editable Button.
$(document).on("click", ".cancelEditableButton", function(event) {
 var oldVal = $(this).siblings(".editField").attr("oldVal");
 $(this).parent().html(oldVal);
 event.stopPropagation();
});


//############  Form Interactions:  ####################################
function getLicensePopup() {
 $("#popupContainer").attr("for", "userLicenseAgreement");
 darkenMask();

 $.get("/user_license_agreement/", function(response) {
  $("#popupContainer_inner").html(response);
  $(".closeButton").hide(); //Remove unnecessary close buttons.
 });

 $("#popupGlobal").fadeIn(300);
}

$(document).on("submit", ".downloadForm", function(event) {
 //Variable Setup
 var form = $(this);
 model=$(this).find("input[name=model]").val();

 loadPageFilters();
 //If filters are enabled, send them in the form.
 var filters = "";
 if ($(this).find("input[name=use_filters]").val()=="True"){
  $(this).find("input[name=filters]").val(JSON.stringify(currentQuery));
 }

 showRibbon("Working...", neutralColor, "#popupContainer");

});

$(document).on("submit", ".infoForm", function() {
 var form = $(this); //Keep a reference to the form inside the POST request.
 if ($(".warningIndicator").length){
  showRibbon("Working...", neutralColor, $("#popupContainer_inner"), false);
 }

 var formName = $(form).attr("name");
 var formAction = $(form).attr("action");

 //Send any get-request params that might need to get through.
 var params = window.location.search;
 var paramsList = params.split("?").slice(1);
 for (var i=0; i < paramsList.length; i++){
   var keyValPair = paramsList[i].split("=");
   $(form).append("<input name=\""+keyValPair[0]+
		"\" type=\"hidden\" value=\""+keyValPair[1]+"\" />");
 }

 var formContents = $(form).serialize();

 $.post(formAction, formContents, function(response) {
  //Remove the loading wheel.
  clearLoadingWheels(form)

  //Translate server-responses to actions.
  if (response=="0"){
   showRibbon("Data added!", goodColor, $("#popupContainer_inner"));
   refreshOnMaskFade();
   return false;
  } else if (response=="0_close"){
   //If the form simply should close, do so.
   refreshOnMaskFade();
   $("#mask").trigger("click");
   return false;
  } else if (response=="1") {
   getLicensePopup();
   return false;
  } else if (response=="2") {
   showRibbon("Edit failed!", badColor, $("#popupContainer_inner"));
   $(form).append("<input type=\"submit\" value=\"Save\" class=\"button\"/>");
   return false;
  } else if (response=="3") {
   showRibbon("Info missing!", badColor, $("#popupContainer_inner"));
   $(form).append("<input type=\"submit\" value=\"Save\" class=\"button\"/>");
   return false;
  } else if (response=="4") {
   showRibbon("Invalid data!", badColor, $("#popupContainer_inner"));
   $(form).append("<input type=\"submit\" value=\"Save\" class=\"button\"/>");
   return false;
  } else if (response=="5") {
   showRibbon("Please select a file!", badColor, $("#popupContainer_inner"));
   $(form).append("<input type=\"submit\" value=\"Upload\" class=\"button\"/>");
   return false;
  } else {
   if ($("#popupContainer_inner").is(":visible")){
    //Recreate the popup window with the server response.
    $("#popupContainer_inner").html(response);
    $(".subPopup").draggable();
   } else {
    $(".infoContainer").html(response);

   }
  }


  //Reset the autocomplete if necessary.
  $(".autocomplete_reactant").autocomplete({
   source: CGAbbrevs,
   messages: {
    noResults: "",
    results: function() {}
   }
  });

  reloadIfNeeded();

  //Show the ribbon message if applicable.
  if ($(".successActivator").length) {
   showRibbon("Added!", goodColor, "#popupContainer_inner");
   refreshOnMaskFade();
   $(".successActivator").remove();
   return false;
  }

 });
 return false; //Do not continue or else the form will post again.
});

$(document).on("click", ".button[type=submit]", function() {
 var form = $(this).closest("form");

 //If form is invalid, display the first error and don't remove the submit button.
 var errors = getFormErrors(form);

 if (errors.length==0){
   $(this).hide();
   applyLoadingWheel($(this).parent("form"));
 } else {
  showRibbon(errors[0], badColor, "body");
  return false;
 }
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

function sendSearchQuery(currentQuery, destination) {
 JSONQuery = JSON.stringify({"currentQuery":currentQuery,"page":"1"});
 $.ajax({
  url:destination,
  method:"post",
  data: {"body":JSONQuery},
  traditional: true,
  success: function(response) {
   var infoBars = $(".infoBar") //Rescue any infoBars that may exist.
   $("#mainPanel").html(response)

   $("#mainPanel").prepend(infoBars)
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

//Set the autocomplete box on reactants if needed.
$(document).on("change", ".dropDownMenu[name=field]", function() {
 if ($(this).val()=="reactant"){
  $(".autocomplete_reactant").autocomplete("enable")
  $("#searchValue").addClass("autocomplete_reactant")
 } else {
  $(".autocomplete_reactant").autocomplete("disable")
 }
});

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
      showRibbon("Already queried!", neutralColor,"#sidePanel", true);
      return false //Don't continue if query is already present.
     }
    }
   }
   showRibbon("Searching!", goodColor,"#sidePanel", true);

   currentQuery.push({
    "field":field,
    "match":"exact",
    "value":value,
   });
   sendSearchQuery(currentQuery, $(this).closest("form").attr("action"));

  //If "Fields" search is active.
  } else if ($(".ui-state-active").children().html()=="Fields" && $("#searchValue").val()){
   field = $(".dropDownMenu[name=field] option:selected").val()
   match = $(".dropDownMenu[name=match] option:selected").val()
   value = $("#searchValue").val()
   //Make sure the query was not already searched.
   if (currentQuery){
    for (var i in currentQuery) {
     if (currentQuery[i]["field"] == field && currentQuery[i]["value"] == value) {
      showRibbon("Already queried!", neutralColor,"#sidePanel", true);
      return false //Don't continue if query is already present.
     }
    }
   }
   showRibbon("Searching!", goodColor,"#sidePanel", true);

   currentQuery.push({
    "field":field,
    "match":match,
    "value":value,
   });
   sendSearchQuery(currentQuery, $(this).closest("form").attr("action"));
  } else {
   showRibbon("Nothing entered!", badColor, $("#sidePanel"), true);
  }
 } else {
  showRibbon("No data to filter!", badColor, $("#sidePanel"), true);
 }
});

$(document).on("click", ".search_backButton", function() {
 if (currentQuery.length){
  currentQuery.pop();
  showRibbon("Removing last filter!", goodColor,"#sidePanel", true);
  sendSearchQuery(currentQuery, $(this).closest("form").attr("action"));
 } else {
  showRibbon("No filters present!", badColor, $("#sidePanel"), true);
 }
});

$(document).on("click", ".search_clearButton", function() {
 currentQuery = Array();
 sendSearchQuery(currentQuery, $(this).closest("form").attr("action"));
 showRibbon("Filters emptied", goodColor,"#sidePanel", true);
});

 //################   Autocompleting Search   ####################### //###
window.setReactantAutoComplete = function() {
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
}

//############ Popup Management: #######################################
window.createPopupConfirmation = function(message) {
 $("body").append("<div class=popupConfirmation title=Confirm:>"+
  message+"</div>");
}

function loadPopup(response, kwargs){
  //Default to no options being set.
  if (kwargs===undefined) kwargs = {};

  //If the user needs to be logged in, send them to a prettier login form.
  if (response.indexOf("/user_login/")>0 && !kwargs["popupLogin"]){
    window.location.href = "/login/";
  } else {
    $("#popupContainer_inner").html(response);

    //Allow popup features on any .autocomplete_reactant element.
    if ($("#popupContainer_inner").find(".autocomplete_reactant").length){
      setReactantAutoComplete();
    }
  }
}

function loadSideBar(response, kwargs){
  //If the user needs to be logged in, send them to a prettier login form.
  if (response.indexOf("/user_login/")>0){
    window.location.href = "/login/";
  } else {
    $("#sidePanel_inner").html(response);
    toggleSideContainer();
    $("#tabs").tabs({active: 1});

    //Allow optional popup features.
    if (kwargs["autocomplete"]){
       setReactantAutoComplete();
    }
  }
}

//Activate specified popup:
$(document).on("click", ".popupActivator", function(event) {
 //Display a loading message if the request takes a visible amount of time.
 $("#popupContainer_inner").html("<center>Loading. Please wait.</center>");

 //Load the activator CSS.
 activatorID = $(this).attr("id");

 switch (activatorID) {
  case "transferRecommendation":
   var pid = $(this).closest(".dataGroup").attr("pid");
   $.get("/data_form/", {pid : pid, model: "rec"}, function(response) {
    loadPopup(response);
   })
   break;
  case "leftMenu_addNew":
   //
   var pid = $(this).closest(".dataGroup").attr("pid");
   //Send the request to the server
   $.get("/data_form/", {pid : pid, model:"data"}, function(response) {
    loadPopup(response, {"autocomplete":true});
   });
   break;
  case "leftMenu_upload_data":
   $.post("/upload_prompt/", {model: "Data"}, function(response) {
    loadPopup(response);
   });
   activatorID = "leftMenu_uploadCSV";
   break;
  case "leftMenu_upload_compoundentry":
   $.post("/upload_prompt/", {model: "CompoundEntry"}, function(response) {
    loadPopup(response);
   });
   activatorID = "leftMenu_uploadCSV";
   break;
  case "addReactantGroup":
   var group = $(this).closest(".reactantField").attr("group");
   var pid = $(this).closest(".dataGroup").attr("pid");
   $.get("/add_reactant/", {group: group, pid:pid}, function(response) {
    loadPopup(response, {"autocomplete":true});
   });
   activatorID = "addReactantGroup";
   event.stopPropagation();
   break;
  case "leftMenu_download_data":
   $.post("/download_prompt/", {model: "Data"}, function(response) {
    loadPopup(response);
   });
   activatorID = "leftMenu_downloadCSV";
   break;
  case "leftMenu_download_compoundentry":
   $.post("/download_prompt/", {model: "CompoundEntry"}, function(response) {
    loadPopup(response);
   });
   activatorID = "leftMenu_downloadCSV";
   break;
  case "leftMenu_download_saved":
   $.post("/download_prompt/", {model: "Saved"}, function(response) {
    loadPopup(response);
   });
   activatorID = "leftMenu_downloadCSV";
   break;
  case "leftMenu_download_recs":
   $.post("/download_prompt/", {model: "Recommendation"}, function(response) {
    loadPopup(response);
   });
   activatorID = "leftMenu_downloadCSV";
   break;
  case "searchButton":
   PT_selected = Array();
   $.get("/search/Data", function(response) {
    loadSideBar(response, {"autocomplete":true})
   });
   return false;
   break;
  case "searchButton_recs":
   $.get("/search/Recommendation", function(response) {
    loadSideBar(response, {"autocomplete":true})
   });
   return false;
   break;
  case "searchButton_seed_recs":
   $.get("/search/SeedRecommendation", function(response) {
    loadSideBar(response, {"autocomplete":true})
   });
   return false;
   break;
  case "userLogin":
   $.get("/user_login/", function(response) {
    loadPopup(response, {"popupLogin":true});
   });
   break;
  case "changePassword":
   $.get("/change_password/", function(response) {
    loadPopup(response);
   });
   break;
  case "registrationPrompt":
   $.get("/registration_prompt/", function(response) {
    loadPopup(response);
   });
   break;
  case "userUpdate":
   $.get("/user_update/", function(response) {
    loadPopup(response);
   });
   break;
  case "userRegistration":
   $.get("/user_registration/", function(response) {
    loadPopup(response);
   });
   break;
  case "labRegistration":
   $.get("/lab_registration/", function(response) {
    loadPopup(response);
   });
   break;
  case "CompoundGuide":
   $.get("/compound_guide/", function(response) {
    loadPopup(response);
    $(".CG_saveButton").remove()
    //Erase the current JSON object and selected CG entries.
    CGEntries = undefined;
    CGAbbrevs = undefined;
    CGSelected = Array();
   });
   break;
  case "CompoundGuide_form":
   $.get("/compound_guide_form/", function(response) {
    loadPopup(response);
    $(".CG_saveButton").remove()
   });
   break;
  default:
   return false; //If a popup is not recognized, don't load anything.
 }
 $("#popupContainer").attr("for", activatorID);
 $("#popupGlobal").fadeIn(300);
});

// Cancel any masks or popups that are covering the screen.
$(document).on("click", ".refreshButton", function() {
 window.location.reload(true);
});

// Grab the rest of the compound_entry to avoid looking up on ChemSpider.
$(document).on("click", ".customCompoundButton", function() {
 var form = $(this).closest("form");

 //Remove any search-related elements from the form.
 $(".customCompoundButton").remove();
 $(".checkCompoundButton").remove();
 $(".compoundCheckResults").remove();
 clearLoadingWheels(form)

 //Show the hidden fields and mark the form as "custom."
 $("input[name=custom]").val(true);
 $(".hiddenField").removeClass("hiddenField");

 if ($("input[type=submit]").length==0){
  $(form).append("<input name=\"submit\" type=\"submit\" value=\"Save\" class=\"button\"/>");
 }

});

// Check a compound before adding it.
$(document).on("click", ".checkCompoundButton", function() {
 //Remove any search-related elements from the form.
 $("input[type=submit]").remove();
 clearLoadingWheels($(this).parent());

 //Variable Setup
 var form = $(this).closest("form");
 var CAS_ID = $(form).find("input[name=CAS_ID]").val();
 var compound = $(form).find("input[name=compound]").val();

 if (CAS_ID || compound) {
  showRibbon("Working...", neutralColor, "#popupContainer", false);
  data = {
    CAS_ID : CAS_ID,
    compound: compound,
   }
  $.get("/check_compound/", data, function(response) {
   $(form).find("input[name=submit]").remove();
   $(".ribbonMessage").remove();
   if (response=="1"){
    showRibbon("Compound not found!", badColor, "#popupContainer_inner");
   } else {
    newResults = "<div>Common Name: "+response["commonName"]+"</div>"
      + "<div>Molecular Formula: "+response["mf"]+"</div>"
      + "<div>Molecular Weight: " + response["mv"]+" g/mol</div>"
      + "<img src=\""+response["imageurl"]+"\" class=\"CG_image_container\"/>"
      + "<div>Is this the compound you want to add?</div>";
    $(".compoundCheckResults").html(formatText(newResults));
    $(form).append("<input name=\"submit\" type=\"submit\" value=\"Save\" class=\"button\"/>");

   }
  });

 } else {
  showRibbon("Nothing to check!", badColor, "#popupContainer_inner");
 }
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

//############ Forced Editable Changes #################################

$(document).on("change", ".editable_assignedUser", function() {
 var JSONQuery = {
  rec_pid:$(this).closest(".dataGroup").attr("pid"),
  user_pid:$(this).children("option:selected").attr("val")
 }

 $.post("/assign_user/", JSONQuery, function(response) {
  if (response==0) {
   showRibbon("Edit successful!", goodColor, "#mainPanel");
  } else {
   showRibbon("Edit failed!", badColor, "#mainPanel");
  }
 });
});

//############ Editable Data: ##############################################
//When the user clicks on an autocomplete option, trigger the "keyup" event.
$(document).on("click", ".ui-menu-item", function() {
 if ($(":focus").attr("class").indexOf("editField") != -1) {
  $(".editText:focus").keyup();
 }
});

//Make edit text fields auto-size and validate while typing.
$(document).on("keyup", ".editText", function(e) {
  // Keyboard shortcuts.
  if (e.keyCode==27) { // The `esc` key
    $(this).siblings(".cancelEditableButton").trigger("click");
  } else if (e.keyCode==13) { // The `enter` key
    $(this).siblings(".editConfirm").trigger("click");
  }

  adaptSize($(this));

  if (parseInt($(this).closest(".reactantField").attr("group"))>2) {
    var required = false;
  } else { var required = true; }

  if (!quickValidate(($(this).parent().attr("class").split(" ")[1]).substr(5),
    $(this).val(), required)) {
    $(this).addClass("badData");
  } else {
    $(this).removeClass("badData");
  }
});

function cancelEditables() {
 $(".editField").each(function() {
  var editParent = $(this).parent(".editable");
  var oldVal = $(this).attr("oldVal");
  $(editParent).html(oldVal);
 });
}

//Initiate edit session.
$(document).on("click", ".editable", function() {
 //Don't allow editables of empties in a reactantField
 if ($(this).closest(".reactantField").length && $(this).is(":empty")){
  return false;
 } else if ($(this).closest(".recommendation").length && $(this).attr("class").indexOf("type_notes")<0) {
  return false;
 }

 //Close any other editables.
 if ($(this).find(".editField").length != 0) {
  return false;
 } else {
  cancelEditables();
 }

 var pidToChange = $(this).closest(".dataGroup").attr("pid");
 $(this).css("opacity",1);

 if ($(this).children(".editConfirm").length == 0 ) {
  var oldVal = String($(this).html());
  var editAs = $(this).attr("editAs");

  if (editAs == "select") {
   var options = getOptions($(this).attr("class").split(" ")[1].substr(5));
   var newInnards = "<select class=\"editField editMenu dropDownMenu\""
    + "pidToChange=\"" + pidToChange + "\" "
    +"oldVal=\""+oldVal+"\">";
   for (var i in options) {
    var choice = options[i];
    if (oldVal == choice) {
     newInnards += "<option value=\""+choice+"\"selected>"+choice+"</option>";
    } else {
     newInnards += "<option value=\""+choice+"\">"+choice+"</option>";
    }
   }
   newInnards += "</select> <input class=\"editConfirm\" type=\"button\" value=\"&#10003;\" />"
   $(this).html(newInnards);
   $(this).children(".editMenu").focus();
   $(this).attr("title","");
  } else {
   var cancelButton = ("<div title=\"Cancel this edit.\""
    +"class=\"genericButton cancelEditableButton\">x</div>");
   var inputFields = ("<input class=\"editField editText\" type=\"" + editAs + "\" "
    + "pidToChange=\"" + pidToChange + "\" "
    + "oldVal=\""+ oldVal + "\" "
    + "value=\""+ oldVal + "\" />"
    + cancelButton
    + "<input class=\"editConfirm\" type=\"button\" value=\"&#10003;\" />"
    );

   $(this).html(inputFields);
   adaptSize($(this).children(".editText"));


   if ($(this).attr("class").indexOf("type_reactant") >= 0) {
    //Get the auto-complete options form the CG guide.
    if (CGAbbrevs == undefined) {
     CGAbbrevs = Array();
     for (var key in CGEntries) {
      CGAbbrevs.push(key);
     }
    }
    $(this).children(".editText").autocomplete({
     source: CGAbbrevs,
     messages: {
      noResults: "",
      results: function() {}
     }
    });
   }
   $(this).children(".editText").focus();
   $(this).attr("title","");
  }
 }
 return false; //Don't continue on to select the data.
});
/* //Needs to check reactants as well. Perhaps just do quick request?
function checkForUnknown(elem) {
 var is_valid = true;
 alert("1");
 $(elem).find(".boolField").each(function() {
  if (is_valid && $(this).html()=="?") {
   $(elem).addClass("badData");
   alert("2");
   return false;
  }
 });
 alert("3");
 $(elem).removeClass("badData");
}*/

//Confirm edit with button press.
$(document).on("click", ".editConfirm", function() {
 var editFieldSibling = $(this).siblings(".editField")
 var editParent = $(this).parent();
 var editRow = $(editParent).parents("tr");
 //Find the general fieldChanged (eg, quantity vs. quantity_1)
 var fieldChanged = $(this).closest(".editable").attr("class").split(" ")[1];
 var numChanged = $(this).closest(".reactantField").attr("group");
 var newValue = $(editFieldSibling).val().trim();
 var oldValue = $(editFieldSibling).attr("oldVal");

 var validData = false;

 try {
  if ($(this).closest(".dataGroup").attr("class").indexOf("recommendation")>=0) {
   var model="rec";
  } else {
   var model="data";
  }
 } catch(err) {}

 if (editFieldSibling.attr("class").split(" ").indexOf("editText") != -1) { //Edit by Text
  //Check if the data is required (ie: if it pertains to reactant 3-5).
  if (parseInt(numChanged)>2) {
   var required = false;
  } else { var required = true; }
  if ($(this).siblings(".editText").attr("class").indexOf("badData") < 0
   && quickValidate(fieldChanged.substr(5), newValue, required)) {
   validData = true;
  } else {
   showRibbon("Invalid!", badColor, "#dataContainer");
  }
 } else { //Edit by Menu
  //Since the only data choices are those which are supplied...
  validData = true;
 }

 if (validData) {
  if (newValue != oldValue) {
   //Find the specific fieldChanged.
   if ($.isNumeric(numChanged)) {
     fieldChanged = fieldChanged.substr(5) + "_" + numChanged;
    } else {
     fieldChanged = fieldChanged.substr(5);
    }

    //Send edits for Compound Guide
    if ($("#CG_display").length) {
     var pidToChange = $(this).closest(".CGRow").attr("pid");
     var compound = $(this).closest(".CGRow").find(".CG_compound div").html();
     showRibbon("Working...", neutralColor, "#popupContainer", false);
     var editLog = {
      field : fieldChanged,
      newVal : newValue,
      oldVal : oldValue,
      pid : pidToChange,
      type : "edit"
     }
     $.post("/edit_CG_entry/", editLog, function(response) {
       if (response!="0"){
        $(editParent).html(oldValue);
        showRibbon(response, badColor, "#popupContainer");
       } else {
         refreshOnMaskFade();
         showRibbon("Entry changed!", goodColor, "#popupContainer");
       }
     });

    } else {
     //Send edits for the Data/Rec View
     var editLog = {
      pid : $(editFieldSibling).attr("pidToChange"),
      field : fieldChanged,
      newValue : newValue,
     };

    //Get the appropriate URL.
    if (model=="rec"){
     var url = "/change_Recommendation/";
    } else {
     var url = "/change_Data/";
    }

    $.post(url, editLog, function(response) {
     //The data should now be up to date:
     if (response != 0) {
      $(editParent).html(oldValue);
      showRibbon(response, badColor, "body");
     } else {
      showRibbon("Entry changed!", goodColor, "#mainPanel");
     }
    });
   }

   //Immediately change the visual for the user while waiting for a response.
   if (typeof compound !== 'undefined') {
    if (fieldChanged=="abbrev" && newValue==compound) {
     $(editParent).empty();
    } else {
     $(editParent).html(newValue);
    }
   } else {
    $(editParent).html(newValue);
    //checkForUnknown($(this).closest(".dataGroup"));
   }

  } else {
   //Revert to old value if new value is unchanged.
   $(editParent).html(oldValue);
  }
  restyleData();
 }
 return false; //Don't re-edit the data (since ".editable" was clicked again).
});

//################ CG Editing: #########################################
//Delete CG data button (but requires a "save" confirmation).
$(document).on("click", ".CG_deleteButton", function() {
 //Add the compound guide entry index to the selected data list.
 var editParent = $(this).closest("tr");

 //Send data to identify the entry to be deleted.
 CGSelected.push($(this).closest(".CGRow").attr("pid"));

 //Display a CG save button if one does not exist.
 $("#popupContainer").append("<div class=\"CG_saveButton genericButton\">Save</div>");
 $("#compoundGuideForm").html("Please save before continuing.");
 //Clear editing abilities and revert any edits-to-be-made.
 $(".editable").removeClass("editable");
 $(".editField").each(function() {
  var oldVal = $(this).attr("oldVal");
  $(this).parent().html(oldVal);
  });

 //Remove the data from the CG visual.
 $(this).closest("tr").remove();
});

$(document).on("click", ".CG_saveButton", function() {
 showRibbon("Working...", neutralColor,"#popupContainer_inner", false);
 //Send the selected CG entry indexes to the server to be deleted.
 var JSONQuery = {type:"del", pids : CGSelected};
 $.post("/edit_CG_entry/", JSONQuery, function(response) {
  if (response==0) {
   //Show a newly updated screen.
   CGSelected = Array();
   window.location.reload(true);
  } else {
   showRibbon("Edit failed!", badColor,"#popupContainer_inner");
  }
 });
});


//############   Upload   ########################################
//Change visible file name when applicable.
$(document).on("change", "#uploadCSV_hiddenInput", function() {
 if ($("#uploadCSV_hiddenInput").val()) {
  $("#uploadCSV_display").children("div").html($("#uploadCSV_hiddenInput").val().split("\\").pop());
 } else {
  $("#uploadCSV_display").children("div").html("None Selected");
 }
});

//############ Post-load Config ########################################
reloadIfNeeded();
restyleData();

//######################################################################
});
