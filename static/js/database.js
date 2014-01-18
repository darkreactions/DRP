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

//###############  Change Data   #########################################
// # # # # # Delete Reactant Button
$(document).on("mouseover", ".reactantField", function() {
 if ($(".reactantButton").length==0 ){
  if ($(this).children(".type_reactant").is(":empty")) {
   $(this).append("<div id=\"addReactantGroup\" class=\"reactantAddButton reactantButton popupActivator genericButton\"" + 
    "title=\"Add this reactant.\">+</div>"); 
  } else {
   $(this).append("<div class=\"reactantRemoveButton reactantButton genericButton\"" + 
    "title=\"Delete this reactant.\">x</div>"); 
  }
 }
});
 
$(document).on("mouseleave", ".reactantField", function() {
 $(".reactantButton").remove();
})

$(document).on("click", ".reactantRemoveButton", function(event) {
 var reactantField = $(this).closest(".reactantField");
 var pid = $(this).closest(".dataGroup").attr("pid");
 var group = $(reactantField).attr("group");

 $(".popupConfirmation").remove();
 createPopupConfirmation("Erase these reactant fields?");
 $(".popupConfirmation").dialog({
  resizable: false,
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
    showRibbon("Working...", neutralColor, "body");
    $.post("/delete_reactant/", {group:group, pid:pid}, function(response) {
     if (response==0){ 
      showRibbon("Data deleted!", goodColor, "body");
      $(reactantField).children(".dataField").empty();
     } else {
      showRibbon(response, badColor, "body");
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

 event.stopPropagation();
})


//Delete Button
//Ask for confirmations for deletions.
$("#leftMenu_delete").click(function() {
 if (selectedData.length) {
  //Delete any extra popup confirmations that exist.
  $(".popupConfirmation").remove();
  createPopupConfirmation("Really delete the " + selectedData.length + " selected data?");
  $(".popupConfirmation").dialog({
   resizable: false,
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

