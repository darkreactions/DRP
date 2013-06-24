$(document).ready(function() {
//############ Variable Setup: #########################################
var selectedData = Array();
var changesMade = {
	del:[],
	edit:[],
	dupl:[],
	add:[],
	};
	
//############ Post-load Setup: #########################################
restyleData() //Style the data according to what it contains.

//############ Dependent Functions: ####################################
//Sort from greatest to least.
function sortNumbers(smallNum,bigNum) {
	return smallNum - bigNum;
}

//Refresh screen.
function refreshScreen(rememberPage) {
	rememberPage = rememberPage !== undefined ? rememberPage : true
	if (rememberPage) { 
		setPageCookie(); 
	}
	window.setTimeout('location.reload()', 300);
}

//Set the current_page cookie.
function setPageCookie(page) {
	page = page !== undefined ? page : $("#pagesCurrent").html().trim();
	$.cookie("current_page", page,
		{expires: 7}); //Set cookie to expire after one week.
}

//Refresh the data container classes. (NOTE: Does not perform server request for new data.)
function restyleData() {
	//Fade out units if the amount is also faded out.
	$(".type_notes").each(function() {
		if ($(this).is(":empty")) {
			$(this).css({"opacity":"0.3"});
		}
	});
	//Keep selected data highlighted even if on page changes.
	$(".dataIndex").each(function() {
		dataID = Number($(this).html().trim())
		if (selectedData.indexOf(dataID) != -1) {
			$(this).parent().addClass("dataSelected");
		}
	});
}

function showRibbon(message, color, location, timeout) {
	//Default location is the dataContainer.
	location = location !== undefined ? location : "#dataContainer"
	//Assume that the ribbon should time out.
	timeout = timeout !== undefined ? timeout : true
	
	//Remove any extra 
	$(".ribbonMessage").remove();
	
	$(location).append(
		"<div class=ribbonMessage style=\"background-color:" + color + 
			";\">" + message + "</div>"
	);
	if (timeout) {
		setTimeout(function() {
			$(location + " .ribbonMessage").fadeOut(1000);
		},500);
	}
}

function createPopupConfirmation(message) {
	$("body").append("<div class=popupConfirmation title=Confirm:>"+
		message+"</div>");
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

//############ Server Transactions: #########################################
//Record a change if the data is valid.
function changeElement(element) { //"element" is the temporary input field.
	var newValue = element.val();
	var originalValue = element.attr("originalText");
	//Remove "type_" from the field name.
	var fieldChanged = element.parent().attr("class").split(' ')[1].substr(5);
	
	//Add the reactant number to the class (eg, 1-5) if applicable.
	if (fieldChanged == "reactant" || fieldChanged == "quantity" || fieldChanged == "unit") {	
		fieldChanged += "_"+ element.parent().attr("class").split(' ')[2];
	}
	
	//If the new value isn't valid, don't ask the server to verify.
	if (clientValidate(fieldChanged, newValue)) {
		var indexChanged = parseInt(element.parent().parent().siblings(
			".dataIndex").html());
		changesMade.edit.push([indexChanged, fieldChanged, newValue]);
	} else {
		showRibbon("Invalid edit!", "red");
		newValue = originalValue; //Revert back to old entry if invalid
	}
	element.parent().html(newValue);
	
	//Only send a POST to the server if a change was made.
	if (originalValue != newValue) {
		submitChanges(); //###Perhaps update every few seconds?
	}	
}

function submitChanges() {
	if ((changesMade.del.length > 10) || (changesMade.dupl.length > 10)) {
		showRibbon("Working. This may take a moment.", "orange","body", false);
	}

	JSONArray = JSON.stringify(changesMade);
	$.post("/data_update/", JSONArray, function() {
		//The data should now be up to date:
		changesMade = {
			del:[],
			edit:[],
			add:[],
			dupl:[],
			}; 
		refreshScreen(true);
	});
}
	
//############ Popup Management: #######################################
//Activate specified popup:
$(document).on("click", ".popupActivator", function() {
	//Display a loading message if the request takes a visible amount of time.
	$("#popupContainer_inner").html("<center>Loading. Please wait.</center>");
	
	//Load the activator CSS.
	activatorID = $(this).attr("id");
	$("#popupContainer").attr("for", activatorID);
	
	switch ($(this).attr("id")) {
		case "dbMenu_addNew":
			$.get("/data_form/", function(response) {
				$("#popupContainer_inner").html(response);
			});
			break;
		case "dbMenu_uploadCSV":
			$.get("/upload_CSV/", function(response) {
				$("#popupContainer_inner").html(response);
			});
			break;
		case "userLogin":
			$.get("/user_login/", function(response) {
				$("#popupContainer_inner").html(response);
			});
			break;
		case "registrationPrompt": //###REPLACE WITH PRELOADED DROP-DOWN MENU?
			$.get("/registration_prompt/", function(response) {
				$("#popupContainer_inner").html(response);
			});
			break;
		case "userRegistration":
			$.get("/user_registration/", function(response) {
				$("#popupContainer_inner").html(response);
			});
			break;
		case "labRegistration":
			$.get("/user_registration/", function(response) {
				$("#labRegistration").html(response);
				
			});
			break;
		case "dbMenu_downloadCSV"://###Replace with cute form?
			$("#popupContainer_inner").html("CSV downloading!");
			//Download the saved data as a CSV.
			location.replace("/download_CSV/");
			break;
		default: 
			return false; //If a popup is not recognized, don't load anything.
	}
	$("#popupGlobal").fadeIn("fast");
	
});

// Fade the popup when the mask is clicked.
$(document).on("click", "#mask", function() {
	$("#popupGlobal").fadeOut("fast");
	
	//Reload the screen if requested.
	if ($(".reloadActivator").length) {
		refreshScreen();
		$(".reloadActivator").remove();
	};
});

//Close popups on close-button click.
$(document).on("click", ".closeButton", function() {
	//Close the global container if the main container is closed.
	if ($(this).parent().attr("id")=="popupContainer") {
		$(this).parent().parent().fadeOut("fast"); 
	} else {
		$(this).parent().fadeOut("fast");
	}
});

//Show the CSV title to be displayed on the visible element.
$(document).on("change", "#uploadCSV_hiddenInput", function() {
	if (!$("#uploadCSV_hiddenInput").val()) {
		$("#uploadCSV_display").children("div").html("None Selected");
	} else {
		$("#uploadCSV_display").children("div").html($(this).val().split('\\').pop());
	}
});

//############ Form Interactions: ######################################

$(document).on("submit", ".uploadForm", function() { //###Ugly...
		//Show the ribbon message if applicable.
		//////if ($(".successActivator").length) {
			showRibbon("Working. This may take a moment.", "orange", "body", false);
			//////$(".successActivator").remove();
		//////}	
});

$(document).on("submit", ".infoForm", function() {
	var form = $(this); //Keep a reference to the form inside the POST request.
	$.post($(form).attr("action"), $(form).serialize(), function(response) {
		//Recreate the popup window with the server response.
		$("#popupContainer_inner").html(response);
		$(".subPopup").draggable();
		
		//Show the ribbon message if applicable.
		if ($(".successActivator").length) {
			showRibbon("Data added!", "green", "#popupContainer_inner");
			$(".successActivator").remove();
			return false;
		}
		
		//Reload the page if applicable.
		if ($(".reloadActivator").length) {
			refreshScreen(false); //Do not remember the log-in page.
		}
		
	});
	
	return false; //Do not continue or else the form will post again. 
});



//////############ Data Selection: #########################################
//Click Group to Select:
$(document).on("click", ".dataGroup", function() {
	$(this).toggleClass("dataSelected");
	
	//Get Index of Selected Data by Number
	var dataID = parseInt($(this).children(".dataIndex").html());
	
	var indexOfData = selectedData.indexOf(dataID);
	if (indexOfData >= 0) { //Data is already selected (thus, deselect)
		selectedData.splice(indexOfData,1);	
	} else { //Select data
		selectedData.push(dataID);
		selectedData.sort(sortNumbers); //###Necessary? Probably
	}
});

//Select Page (Button)
$("#dbMenu_selectPage").click(function() { 
	var selectedOnPage = [];
	
	//Select data
	$(".dataIndex").each(function() {
		var dataID = parseInt($(this).html());
		if (selectedData.indexOf(dataID) < 0){ //If the item is not in the list yet.
			$(this).parent().addClass("dataSelected");
			selectedData.push(dataID);
		} else {
			//Count the element as selected.
			selectedOnPage.push(dataID);
		} 
	}); 
	
	//Deselect data (if all elements were already selected; ie, "toggle")
	if (selectedOnPage.length == $(".dataGroup").length){
		for (i in selectedOnPage) {
			var indexOfData = selectedData.indexOf(selectedOnPage[i]);
			selectedData.splice(indexOfData,1);	
		}
		$(".dataGroup").removeClass("dataSelected");
	}
	selectedData.sort(sortNumbers); //Sort data once selection is complete.
});

//Select All (Button)
$("#dbMenu_selectAll").click(function() { 
	//If all data is selected, deselect everything.
	var totalDataSize = $("#pagesTotal").attr("total_data_size");
	if (selectedData.length == totalDataSize) {
		selectedData = [];
		$(".dataGroup").each(function() {
			$(this).removeClass("dataSelected");
		});
	} else {
		selectedData = [] //Clear selection to avoid duplicates.
		for(i = 1; i <= totalDataSize; i++) { 
			selectedData.push(i);
		}
		//Data will already be sorted from least to greatest after loop. 
		$(".dataGroup").addClass("dataSelected");
	}
});

//View Full Datum (Button)
$(document).on("click", ".expandButton", function() {
	//Send the 0-based dataIndex to the server and replace the dataEntry.
	var indexRequested = parseInt($(this).siblings(".dataIndex").html())-1;
	var moreButton = $(this);
	$(moreButton).siblings(".dataEntry").html("Loading. Please Wait.");
	$.post("/get_full_datum/", {"indexRequested":indexRequested}, 
		function(response) {
			$(moreButton).siblings(".dataEntry").html(response);
			$(moreButton).fadeOut("slow");
	});
	
	return false; //Do not continue so the data is not selected.
});

//############### Change Data: #########################################
//Duplicate Button
$("#dbMenu_duplicate").click(function() {  
	if (selectedData.length) { 
		for (var i in selectedData) {
			changesMade.dupl.push(selectedData[i]);
		}
		
		showRibbon("Selection duplicated!", "green");
		
		//Upload updated data
		submitChanges();
	}
});

//Delete Button
//Ask for confirmations for deletions.
$("#dbMenu_delete").click(function() { 
	if (selectedData.length) { 
		//Delete any extra popup confirmations that exist.
		$(".popupConfirmation").remove();
		createPopupConfirmation("Really delete the " + selectedData.length + " selected data?");
		$(".popupConfirmation").dialog({
			resizable: false,
			closeOnEscape: false, // ###Temporary fix?
			height: 150,
			modal: true,
			buttons: {
				"Delete Selection": function() {
					for (var i in selectedData) {
						changesMade.del.push(selectedData[i]);
						//Immediately delete the data from the client's view ("#g_x" represents group ID).
						$("#g_"+String(selectedData[i])).remove()
					}
					selectedData = [];
					
					showRibbon("Selection deleted!", "green");
					
					//Upload updated data
					submitChanges();
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


//############ Edit Data: ###########################@##################
//Initiate edit session.
$(document).on("click", ".editable", function() {
	$(".editable").css("opacity",1);
	
	if ($(this).children(".editConfirm").length == 0 ) {
		var oldVal = String($(this).html());
		var editAs = $(this).attr("editAs");
		$(this).html("<input class=\"editField\" type=\"" + editAs
			+ "\" title=\"Enter the new value.\" oldVal=\""+ oldVal
			+ "\" value=\"" + oldVal + "\" />"
			+ "<input class=\"editConfirm\" type=\"button\" value=\"OK\" />"
			);
		adaptSize($(this).children(".editField"));
		$(this).children(".editField").focus();
	}
	return false; //Don't continue on to select the data.
});

//Confirm edit with button press.
$(document).on("click", ".editConfirm", function() {
	//Validate data and revert to old value if new value is invalid.
	var editFieldSibling = $(this).siblings(".editField")
	var newValue = $(editFieldSibling).val();
	var oldValue = $(editFieldSibling).attr("oldVal");
	if ((true) && (newValue != oldValue)) {
		//Find the fieldChanged.
		var fieldChanged = $(this).closest(".editable").attr("class").split(' ');
		if ($.isNumeric(fieldChanged[2])) { 
			fieldChanged = fieldChanged[1].substr(5) + "_" + fieldChanged[2];
		} else {
			fieldChanged = fieldChanged[1].substr(5);
		}
		
		//Submit the new value to the server. 
		changesMade.edit.push([ //[indexChanged, fieldChanged, newValue]
			$(this).parent().parent().siblings(".dataIndex").html().trim(),
			fieldChanged,
			newValue
			]);
		
		//Immediately change the visual for the user.
		$(this).parent().html(newValue);
	} else {
		$(this).parent().html(oldValue);
	}
	restyleData();
	//submitChanges();
	return false; //Don't re-edit the data (since ".editable" was clicked again).
});

//Make edit text fields auto-size while typing.
$(document).on("keyup", ".editField", function() {
	adaptSize($(this));
});

//############ User Authentication: ####################################
$("#userLogOut").click( function() {
	//Send the log-out signal.
	$.get("/user_logout/", function() {});
	
	//Reload the screen to verify log-out.
	refreshScreen()
});

//############## Change Pages: #########################################
$(document).on("click", "#pagesInputButton", function() {
	var pageLink = $("#pagesInputText").val().trim();
	
	if (pageLink == $("#pageLinkCurrent").html().trim()) {
		showRibbon("Already here!", "green");
		return false //Don't do anything else if page already is active.
	} else if ($.isNumeric(pageLink)
		&& parseInt(pageLink) <= parseInt($("#pagesTotal").html().trim())
		&& parseInt(pageLink) > 0) {
		//Create a CSRF Token for Django authorization.
		var csrftoken = $.cookie("csrftoken");
		
		//Request the information about a page.
		pageDestination = "/data_transmit/" + pageLink;
		$.get(pageDestination, function(dataResponse) {
			var newDataContainer = dataResponse;
			$("#dataContainer").html(newDataContainer);
			setPageCookie(pageLink)
			restyleData();
		});
	} else {
		showRibbon("Page does not exist", "red");
	}
});

$(document).on("click", ".pageLink", function() {
	//If a valid page link has been clicked.
	var pageLink = $(this).html().trim();
	//If the link is a number (not an ellipsis) and not the current page:
	if ($.isNumeric(pageLink) && pageLink != $("#pageLinkCurrent").html().trim()) {
		//Request the information about a page.
		pageDestination = "/data_transmit/" + pageLink;
		$.get(pageDestination, function(dataResponse) {
			var newDataContainer = dataResponse;
			$("#dataContainer").html(newDataContainer);
			setPageCookie(pageLink)
			restyleData();
		});
	}
});


//######################################################################
});
