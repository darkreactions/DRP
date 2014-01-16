$(document).ready(function() {
//############ Variable Setup: #########################################
var CGSelected = Array();

//############ Dependent Functions: ####################################
//Set the current_page cookie.
function setPageCookie(page) {
 page = page !== undefined ? page : $("#pagesCurrent").html().trim();
 $.cookie("current_page", page,
  {expires: 7}); //Set cookie to expire after one week.
}

function adaptSize(element) {
 var numChars = $(element).val().length;
 var maxWidth = $(element).css("max-width") !== "none" ? parseInt($(element).css("max-width")) : 100;
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

//############ CSV Management: #######################################
//Change visible file name when applicable.
$(document).on("change", "#uploadCSV_hiddenInput", function() {
 if ($("#uploadCSV_hiddenInput").val()) {
  $("#uploadCSV_display").children("div").html($("#uploadCSV_hiddenInput").val().split("\\").pop());
 } else {
  $("#uploadCSV_display").children("div").html("None Selected");
 }
});

//Hide the Download Form's dataFilter selection menu if not applicable.
$(document).on("click", "#radio_reactionData", function() {
 $(".dropDownMenu").html(
 "<option value=\"all\">All</option>"+
 "<option value=\"good\">Valid</option>"+
 "<option value=\"bad\">Invalid</option>");
});

$(document).on("click", "#radio_compoundGuide", function() {
 $(".dropDownMenu").html(
 "<option value=\"simple\">Simple</option>"+
 "<option value=\"complex\">Complex</option>");
});


//############ Form Interactions: ######################################
$(document).on("submit", ".uploadForm", function() { //###Ugly...
 showRibbon("Working! This may take a moment.", "#FFC87C", "body", false);
});

//################ CG Editing: #########################################
//Delete CG data button (but requires a "save" confirmation).
$(document).on("click", ".CG_deleteButton", function() {
 //Add the compound guide entry index to the selected data list.
 var editParent = $(this).closest("tr");

 //Send data to identify the entry to be deleted.
 CGSelected.push({
  "compound": $(editParent).find(".type_compound").html(),
  "abbrev": $(editParent).find(".type_abbrev").html()
  });

 //Display a CG save button if one does not exist.
 $("#popupContainer").append("<div class=\"CG_saveButton\">Save</div>");
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
 showRibbon("Saving!", "99FF5E","#popupContainer_inner", false);
 //Sort the list prior to sending.
 CGSelected.sort(sortNumbersReverse);
 //Send the selected CG entry indexes to the server to be deleted.
 JSONArray = JSON.stringify({"type":"del", "data":CGSelected});
 $.post("/edit_CG_entry/", JSONArray, function() {
  //Show a newly updated screen.
  CGSelected = Array();
  window.location.reload(true);
 });
});

//############ Data Selection: #########################################
//Click Group to Select:
$(document).on("click", ".dataGroup", function() {
 $(this).toggleClass("dataSelected");

 //Get Index of Selected Data by Number
 var dataID = $(this).find(".dataEntry>.type_ref").html().trim();
 var indexOfData = selectedData.indexOf(dataID);
 if (indexOfData >= 0) { //Data is already selected (thus, deselect)
  selectedData.splice(indexOfData,1);
 } else { //Select data
  selectedData.push(dataID);
  selectedData.sort(sortNumbers); //###Necessary? Probably
 }
});

//Select Page (Button)
$("#leftMenu_selectPage").click(function() {
 var selectedOnPage = [];
 //Select data
 $(".type_ref").each(function() {
  var dataRef = $(this).html().trim();
  if (selectedData.indexOf(dataRef) < 0){ //If the item is not in the list yet.
   $(this).closest(".dataGroup").addClass("dataSelected");
   selectedData.push(dataRef);
  } else {
   //Count the element as selected.
   selectedOnPage.push(dataRef);
  }
 });

 //Deselect data (if all elements were already selected; ie, "toggle")
 if (selectedOnPage.length == $(".dataEntry").length){
  for (var i in selectedOnPage) {
   var indexOfData = selectedData.indexOf(selectedOnPage[i]);
   selectedData.splice(indexOfData,1);
  }
  $(".dataGroup").removeClass("dataSelected");
 }
});

//############### Change Data: #########################################

//Delete Button
//Ask for confirmations for deletions.
$("#leftMenu_delete").click(function() {
 if (selectedData.length) {
  //Delete any extra popup confirmations that exist.
  $(".popupConfirmation").remove();
  createPopupConfirmation("Really delete the " + selectedData.length + " selected data?");
  $(".popupConfirmation").dialog({
   resizable: false,
   closeOnEscape: false, // ###Temporary fix?
   height: 150,
   position:{
    my:"center",
    at:"center",
    of:"body",
   },
   modal: true,
   buttons: {
    "Delete Selection": function() {
     //Upload updated data
     JSONArray = JSON.stringify(selectedData);
     $.post("/delete_Data/", JSONArray, function(response) {
      showRibbon(response, "#99FF5E", "body");
      if (response == 0) {
       showRibbon("Working...", "#FFC87C", "body", false);
       window.location.reload(true);
      } else {
       showRibbon(response,"#FF6870", "body");
      }
     });

     $(this).dialog("close");
     $(this).remove();
    },
    "No": function() {
     $(this).dialog("close");
     $(this).remove();
    }
   }
  });
 }
});

//############ Edit Data: ##############################################
function cancelEditables() {
 $(".editField").each(function() {
  var editParent = $(this).parent(".editable");
  var oldVal = $(this).attr("oldVal");
  $(editParent).html(oldVal);
 });
}


//Initiate edit session.
$(document).on("click", ".editable", function() {

 //Close any other editables.
 if ($(this).find(".editField").length != 0) {
  return false;
 } else {
  cancelEditables();
 }

 var refToChange = $(this).parent().children(".type_ref").html().trim();

 $(this).css("opacity",1);

 if ($(this).children(".editConfirm").length == 0 ) {
  var oldVal = String($(this).html());
  var editAs = $(this).attr("editAs");

  if (editAs == "select") {
   var options = getOptions($(this).attr("class").split(" ")[1].substr(5));
   var newInnards = "<select class=\"editField editMenu dropDownMenu\""
    + "refToChange=\"" + refToChange + "\" "
    +"oldVal=\""+oldVal+"\">";
   for (var i in options) {
    var choice = options[i];
    if (oldVal == choice) {
     newInnards += "<option value=\""+choice+"\"selected>"+choice+"</option>";
    } else {
     newInnards += "<option value=\""+choice+"\">"+choice+"</option>";
    }
   }
   newInnards += "</select> <input class=\"editConfirm\" type=\"button\" value=\"OK\" />"
   $(this).html(newInnards);
   $(this).children(".editMenu").focus();
   $(this).attr("title","");
  } else {
   $(this).html("<input class=\"editField editText\" type=\"" + editAs + "\" "
    + "refToChange=\"" + refToChange + "\" "
    + "oldVal=\""+ oldVal + "\" "
    + "value=\""+ oldVal + "\" />"
    + "<input class=\"editConfirm\" type=\"button\" value=\"OK\" />"
    );
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

//Confirm edit with button press.
$(document).on("click", ".editConfirm", function() {
 var editFieldSibling = $(this).siblings(".editField")
 var editParent = $(this).parent();
 var editRow = $(editParent).parents("tr");
 //Find the general fieldChanged (eg, quantity vs. quantity_1)
 var fieldChanged = $(this).closest(".editable").attr("class").split(" ")[1];
 var numChanged = $(this).closest(".editable").attr("group");
 var newValue = $(editFieldSibling).val().trim();
 var oldValue = $(editFieldSibling).attr("oldVal");

 var validData = false;

 if (editFieldSibling.attr("class").split(" ").indexOf("editText") != -1) { //Edit by Text
  //Check if the data is required (ie: if it pertains to reactant 3-5).
  if (parseInt(numChanged)>2) {
   var required = false;
  } else { var required = true; }
  if ($(this).siblings(".editText").attr("class").indexOf("badData") < 0
   && quickValidate(fieldChanged.substr(5), newValue, required)) {
   validData = true;
  } else {
   showRibbon("Invalid!", "#FF6870", "#dataContainer");
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
    if (fieldChanged=="compound"){
     var compound = oldValue;
    } else {
     var compound_container = $(this).parent().parent().siblings(".CG_compound").children()
     var compound = $(compound_container).html();
    }

    showRibbon("Working...", "#FFC87C", "#popupContainer", false);
    $.post("/edit_CG_entry/", JSON.stringify({
      "field" : fieldChanged,
      "newVal" : newValue,
      "oldVal" : oldValue,
      "compound": compound, //Acts as the reference to know which Entry to change.
      "type" : "edit"
     }), function(response) {
      if (response!="0"){
       $(editParent).html(oldValue);
       showRibbon(response, "#FF6870", "#popupContainer");
      } else {
       if (fieldChanged=="compound") { //If the compound itself was changed, search for the new value.
        compound = newValue;
       }
       $.post("/compound_guide_entry/", JSON.stringify({"compound":compound}),
        function(response) {
         $(editRow).html(response);
         showRibbon("Success!", "#99FF5E", "#popupContainer");
         //if ($(editRow).find(".type_abbrev").is(':empty')) {
         // $(editRow).find(".type_abbrev").html(oldValue);
         //}
       });
      }
    });

   } else {
    //Send edits for the Database View
    editLog = {
     "ref":$(editFieldSibling).attr("refToChange"),
     "field":fieldChanged,
     "newValue":newValue,
    };
   }
   //Immediately change the visual for the user.
   $(this).parent().html(newValue);
   JSONArray = JSON.stringify(editLog);
   $.post("/change_Data/", JSONArray, function(response) {
    //The data should now be up to date:
    if (response != 0) {
     $(editParent).html(oldValue);
     showRibbon(response, "#FF6870", "body");
   }
   });
   
  } else {
   //Revert to old value if new value is unchanged.
   $(editParent).html(oldValue);
  }
  restyleData();
 }
 return false; //Don't re-edit the data (since ".editable" was clicked again).
});

//When the user clicks on an autocomplete option, trigger the "keyup" event.
$(document).on("click", ".ui-menu-item", function() {
 if ($(":focus").attr("class").indexOf("editField") != -1) {
  $(".editText:focus").keyup();
 }
});

//Make edit text fields auto-size and validate while typing.
$(document).on("keyup", ".editText", function() {
 adaptSize($(this));
 if (parseInt($(this).parent().attr("group"))>2) {
  var required = false;
 } else { var required = true; }

 if (!quickValidate(($(this).parent().attr("class").split(" ")[1]).substr(5),
  $(this).val(), required)) {
  $(this).addClass("badData");
 } else {
  $(this).removeClass("badData");
 }
});

//############## Edit CG: #########################################
//Click to select entire compound.
$(document).on("click", ".CG_compound", function() {

});

//############## Change Pages: #########################################
function changePageTo(page) {
 showRibbon("Loading!", "#FFC87C", "#mainPanel", false);

 //Send the currentQuery if it exists.
 JSONQuery = JSON.stringify({"currentQuery":currentQuery, "page":page});

 $.post("/data_transmit/", {"body":JSONQuery}, function(response) {
  $("#mainPanel").html(response);
  if ($(".dataGroup").length) {
   setPageCookie(page)
   restyleData();
  $(".ribbonMessage").remove();
  }
 });
//TODO: add this back!
// } else {
//  showRibbon("Page does not exist", "#FF6870", "body");
// }
}

$(document).on("click", "#pagesInputButton", function() {
 var pageNum = $("#pagesInputText").val().trim();

 if (pageNum == $("#pageLinkCurrent").html().trim()) {
  showRibbon("Already here!", "#99FF5E", "#dataContainer");
  return false //Don't do anything else if page already is active.
 } else if ($.isNumeric(pageNum)
  && parseInt(pageNum) <= parseInt($("#pagesTotal").html().trim())
  && parseInt(pageNum) > 0) {
  //Create a CSRF Token for Django authorization.
  var csrftoken = $.cookie("csrftoken");
  changePageTo(pageNum)
  }
});

$(document).on("click", ".pageLink", function() {
 //If a valid page link has been clicked.
 var pageNum = $(this).html().trim();
 //If the link is a number (not an ellipsis) and not the current page:
 if ($.isNumeric(pageNum) && pageNum != $("#pageLinkCurrent").html().trim()) {
  //Request the information about a page.
  changePageTo(pageNum)
 }
});


//######################################################################
});

