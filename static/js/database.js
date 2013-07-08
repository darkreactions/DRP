$(document).ready(function() {
//############ Variable Setup: #########################################
var selectedData = Array();
var changesMade = {
	del:[],
	edit:[],
	dupl:[],
	add:[],
	};
var CGSelected = Array();
var CGEntries;
var CGAbbrevs;
	
//############ Post-load Setup: #########################################
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
restyleData() //Style the data according to what it contains.

//############ Tooltip functionality: ##################################
$(document).on("mouseover", ".type_reactant", function() {
	if ($(this).children("input").length == 0){
		if (CGEntries == undefined) {
			$(this).attr("title", "Loading!")
			$.get("/send_CG_names/", function(response) {
				CGEntries = response;
				var compound = CGEntries[$(this).html()] || "Compound not in guide!"
				$(this).attr("title", compound);
				$(".ui-tooltip").html(compound);
				});
		} else {
			var compound = CGEntries[$(this).html()] || "Compound not in guide!"
			$(this).attr("title", compound)
		}
	}
});

var tooltips_disabled = false;
$(document).on("focusin", "input", function() {
	if (!tooltips_disabled) {
		tooltips_disabled = true
		$(document).tooltip("disable");
	}
});

$(document).on("focusout", "input", function() {
	tooltips_disabled = false
	$(document).tooltip("enable");
});

//############ Dependent Functions: ####################################
//Sort from greatest to least.
function sortNumbers(smallNum, bigNum) {
	return smallNum - bigNum;
}

function sortNumbersReverse(smallNum, bigNum) {
	return bigNum - smallNum;
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

//Get edit-by-menu options from editChoices.json.
function getOptions(field) {
	switch (field) {
		case ("unit"):
			return  editChoices["unitChoices"]
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

//############ Server Transactions: #########################################
function submitChanges(refresh) {
	refresh = refresh !== undefined ? refresh : true 
	
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
		if (refresh) {
			refreshScreen();
		}
	});
}

function getCGEntries() {
	
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
				
				//Get the auto-complete options form the CG guide.
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

//################ CG Editing: #########################################
//Delete CG data button (but requires a "save" confirmation).
$(document).on("click", ".CG_deleteButton", function() {
	//Add the compound guide entry index to the selected data list.
	var CGIndex = (parseInt($(this).attr("id").substr(3))-1);
	CGSelected.push(CGIndex);
	
	//Display a CG save button if one does not exist.
	$("#popupContainer").append("<div class=\"CG_saveButton genericButton\">Save</div>");
	$("#compoundGuideForm").html("Please save before continuing.")
	
	//Remove the data from the CG visual.
	$(this).closest("tr").remove();
});

$(document).on("click", ".CG_saveButton", function() {
	showRibbon("Submitting Changes.", "green","#CG_display", false);
	//Sort the list prior to sending.
	CGSelected.sort(sortNumbersReverse);

	//Send the selected CG entry indexes to the server to be deleted. 
	JSONArray = JSON.stringify(CGSelected);
	$.post("/edit_CG_entry/", JSONArray, function() {
		//Show a newly updated screen.
		CGSelected = Array(); 
		refreshScreen(false);
	});
});

//############ Data Selection: #########################################
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
			restyleData();
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
			position:{
				my:"center",
				at:"center",
				of:"#dataContainer_inner",
			},
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
		
		if (editAs == "select") {
			var options = getOptions($(this).attr("class").split(" ")[1].substr(5));
			var newInnards = "<select class=\"editField editMenu dropDownMenu\""
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
			$(this).html("<input class=\"editField editText\" type=\"" + editAs
				+ "\" oldVal=\""+ oldVal
				+ "\" value=\"" + oldVal + "\" />"
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
	//Find the general fieldChanged (eg, quantity vs. quantity_1)
	var fieldChanged = $(this).closest(".editable").attr("class").split(' ');
	var newValue = $(editFieldSibling).val();
	var oldValue = $(editFieldSibling).attr("oldVal");
	
	var validData = false;
	
	if (editFieldSibling.attr("class").split(" ")[1] == "editText") { //Edit by Text
		if ($(this).siblings(".editText").attr("class").indexOf("badData") < 0 
			&& fullValidate(fieldChanged[1].substr(5), newValue)) {
			validData = true;
		} else {
			showRibbon("Invalid data!", "red");
		}
	} else { //Edit by Menu
		//Since the only data choices are those which are supplied...
		validData = true;
	}
	
	if (validData) {
		if (newValue != oldValue) {
			//Find the specific fieldChanged.
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
			submitChanges(false);
		} else {
			//Revert to old value if new value is unchanged.
			$(this).parent().html(oldValue);
		}
		restyleData();
	}
	return false; //Don't re-edit the data (since ".editable" was clicked again).
});

//Make edit text fields auto-size and validate while typing.
$(document).on("keyup", ".editText", function() {
	adaptSize($(this));
	if (!quickValidate(($(this).parent().attr("class").split(" ")[1]).substr(5),
		$(this).val())) {
		$(this).addClass("badData");
	} else {
		$(this).removeClass("badData");
	}
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
			if ($(".dataGroup").length) {
				setPageCookie(pageLink)
				restyleData();
			}
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
