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

//############ Post-load Setup: #########################################
restyleData() //Style the data according to what it contains.

//############ Reactant Tooltips: ##################################
$(document).on("mouseover", ".type_reactant", function() {
	if ($(this).is(":empty")) {
		$(this).attr("title", "Click to add a new reactant.")
	} else if ($(this).children("input").length == 0 && !(tooltips_disabled)){
		if (CGEntries == undefined) {
			var specificDiv = $(this);
			$(specificDiv).attr("title", "Loading!")
			$.get("/send_CG_names/", function(response) {
				CGEntries = response;
				var compound = CGEntries[$(specificDiv).html()] || "Compound not in guide!"
				$(".ui-tooltip").html(compound);
				});
		} else {
			var compound = CGEntries[$(this).html()] || "Compound not in guide!"
			$(this).attr("title", compound)
		}
	} else {
		$(this).attr("title", "")
	}
});

//############ Dependent Functions: ####################################
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

//############ Server Transactions: #########################################
function submitChanges(refresh) {
	refresh = refresh !== undefined ? refresh : true

	if ((changesMade.del.length > 10) || (changesMade.dupl.length > 10)) {
		showRibbon("Working! This may take a moment.", "#FFC87C", "body", false);
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
			window.location.reload(true);
		}
	});
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
	$("#filterContainer").animate({"opacity": 1}, 100);
});
$(document).on("click", "#radio_compoundGuide", function() {
	$("#filterContainer").animate({"opacity": 0}, 100);
});


//############ Form Interactions: ######################################
$(document).on("submit", ".uploadForm", function() { //###Ugly...
	showRibbon("Working! This may take a moment.", "#FFC87C", "body", false);
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
	showRibbon("Saving!", "99FF5E","#popupContainer_inner", false);
	//Sort the list prior to sending.
	CGSelected.sort(sortNumbersReverse);
	//Send the selected CG entry indexes to the server to be deleted.
	JSONArray = JSON.stringify(CGSelected);
	$.post("/edit_CG_entry/", JSONArray, function() {
		//Show a newly updated screen.
		CGSelected = Array();
		window.location.reload(true);
	});
});

//############ Data Selection: #########################################
//Click Group to Select:
$(document).on("click", ".dataGroup", function() {
	//Don't select data in search containers.
	if ($(this).parent().attr("id") == "searchContainer") {
		return false;
	}

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
$("#leftMenu_selectPage").click(function() {
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
$("#leftMenu_selectAll").click(function() {
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

//Create the "Duplicate This Data" button.
$(document).on("mouseover", ".dataGroup", function() {
	if ($(".duplicateSpecificDataButton").length == 0 ){
		var buttonDiv = "<div id=\"leftMenu_addNew_copy\"class=\"";
		buttonDiv += "duplicateSpecificDataButton popupActivator";
		buttonDiv += " genericButton\" style=\"background-image: url(";
		buttonDiv += STATIC_URL+"/icons/add.png);\" title=\""
		buttonDiv += "Copy this reaction to the data form."
		buttonDiv += "\"></div>";
		$(this).append(buttonDiv);
	}
});

$(document).on("mouseleave", ".dataGroup", function() {
	//Eliminate any duplication buttons on when the mouse isn't on a group.
	$(".duplicateSpecificDataButton").remove();
});

//############### Change Data: #########################################

//Duplicate Button
$("#leftMenu_duplicate").click(function() {
	if (selectedData.length) {
		for (var i in selectedData) {
			changesMade.dupl.push(selectedData[i]);
		}
		showRibbon("Duplicated!", "#99FF5E", "#dataContainer", true);

		//Upload updated data
		submitChanges();
	}
});

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

					showRibbon("Deleted!", "#99FF5E", "#dataContainer");

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
var editExemptions = ["searchResultsContainer"]
$(document).on("click", ".editable", function() {
	//Skip any data that should not be edited. //###
	if (editExemptions.indexOf($(".editable").parent().parent().attr("id")) != -1) {
		return false;
	}

	$(this).css("opacity",1);

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
		if (parseInt($(this).parent().attr("class").split(" ")[2])>2) {
			var required = false;
		} else { var required = true; }

		if ($(this).siblings(".editText").attr("class").indexOf("badData") < 0
			&& quickValidate(fieldChanged[1].substr(5), newValue, required)) {
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
			if ($.isNumeric(fieldChanged[2])) {
				fieldChanged = fieldChanged[1].substr(5) + "_" + fieldChanged[2];
			} else {
				fieldChanged = fieldChanged[1].substr(5);
			}

			//Send edits for Compound Guide
			if ($("#CG_display").length) {
				$.post("/edit_CG_entry/", JSON.stringify({
						"field" : fieldChanged,
						"newVal" : newValue,
						"oldVal" : oldValue,
						"type" : "edit"
					}), function() {
				alert("done");
				});

			} else {
				//Send edits for the Database View
				changesMade.edit.push([ //[indexChanged, fieldChanged, newValue]
					$(this).parent().parent().siblings(".dataIndex").html().trim(),
					fieldChanged,
					newValue
					]);
			}

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

//When the user clicks on an autocomplete option, trigger the "keyup" event.
$(document).on("click", ".ui-menu-item", function() {
	if ($(":focus").attr("class").indexOf("editField") != -1) {
		$(".editText:focus").keyup();
	}
});

//Make edit text fields auto-size and validate while typing.
$(document).on("keyup", ".editText", function() {
	adaptSize($(this));
	if (parseInt($(this).parent().attr("class").split(" ")[2])>2) {
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
$(document).on("click", "#pagesInputButton", function() {
	var pageLink = $("#pagesInputText").val().trim();

	if (pageLink == $("#pageLinkCurrent").html().trim()) {
		showRibbon("Already here!", "#99FF5E", "#dataContainer");
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
		showRibbon("Page does not exist", "#FF6870", "#dataContainer");
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
