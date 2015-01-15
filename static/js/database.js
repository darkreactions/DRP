$(document).ready(function() {
//############ Variable Setup: #########################################
var CGSelected = Array();


//############ Data Selection: #########################################
//Click Group to Select:
$(document).on("click", ".dataGroup", function() {
 var edits = $(this).find(".cancelEditableButton")
 if (edits.length) {
   edits.trigger("click");
   return false;
 }

 $(this).toggleClass("dataSelected");


 //Get Index of Selected Data by Number
 var dataID = $(this).attr("pid");
 var indexOfData = selectedData.indexOf(dataID);
 if (indexOfData >= 0) { //Data is already selected (thus, deselect)
  selectedData.splice(indexOfData,1);
 } else { //Select data
  selectedData.push(dataID);
 }
});

//Select Page (Button)
$("#leftMenu_selectPage").click(function() {
 var selectedOnPage = [];
 //Select data
 $(".dataGroup").each(function() {
  var dataPID = $(this).attr("pid");
  if (selectedData.indexOf(dataPID) < 0){ //If the item is not in the list yet.
   $(this).addClass("dataSelected");
   selectedData.push(dataPID);
  } else {
   //Count the element as selected.
   selectedOnPage.push(dataPID);
  }
 });

 //Deselect data (if all elements were already selected; ie, allow "toggle")
 if (selectedOnPage.length == $(".dataGroup").length){
  selectedData = Array();
  $(".dataGroup").removeClass("dataSelected");
 }
});

//###############  Change Data   #########################################
// # # # # # Delete Reactant Button
$(document).on("mouseover", ".reactantField", function() {
 if ($(".reactantButton").length==0 ){
  if ($(this).children(".type_reactant").is(":empty")) {
   $(this).append("<div id=\"addReactantGroup\" class=\"reactantAddButton reactantButton popupActivator genericButton\"" +
    "title=\"Add a reactant here.\">+</div>");
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

});
