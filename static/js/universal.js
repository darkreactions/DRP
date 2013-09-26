$(document).ready(function() {
//######################################################################
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

//############   Ribbons   #############################################
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
			$(location).children(".ribbonMessage").fadeOut(1000);
		},500);
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


	});
	return false; //Do not continue or else the form will post again.
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

function sendSearchQuery(current_query) {
	$.ajax({
		url:"/search/",
		method:"post",
		data: {"current_query":JSON.stringify(current_query)},
		traditional: true,
		success: function(response) {
			$("#searchResultsOuterContainer").html(response)
			if ($("#searchResultsContainer").html().trim()=="No data found!"){
				current_query = []; //Clear the current query if nothing is found.
			}
		}
	});
}

//Back button tooltip
$(document).on("mouseover", "#search_backButton", function() {
	if (current_query.length){
		//Get the previously used filters.
		var filter_string = "Filters:"
		for (var i in current_query.slice(0,-1)) {
			filter_string += "<br/>"+(parseInt(i)+1)+".) "+make_name_verbose(current_query[i]["field"])+": " + current_query[i]["value"]
		}

		//Add the new filter.
		filter_string += "<br/><div class=\"search_backText\">"+parseInt(current_query.length)+".) "+make_name_verbose(current_query[current_query.length-1]["field"])+": "+current_query[current_query.length-1]["value"]+"</div>"

		$(this).attr("title", filter_string);
	} else {
		$(this).attr("title", "");
	}
});

//Filter button tooltip
$(document).on("mouseover", "#search_filterButton", function() {
	if ($("#searchValue").val()){
		//Get the previously used filters.
		var filter_string = "Filters:"
		for (var i in current_query) {
			filter_string += "<br/>"+(parseInt(i)+1)+".) "+make_name_verbose(current_query[i]["field"])+": " + current_query[i]["value"]
		}

		//Add the new filter.
		filter_string += "<br/><div class=\"search_filterText\">"+(parseInt(current_query.length)+1)+".) "+make_name_verbose($("input[name=field]:checked").val())+": "+$("#searchValue").val()+"</div>"

		$(this).attr("title", filter_string);
	} else {
		$(this).attr("title", "");
	}
});

$(document).on("click", "#search_filterButton", function() {
	if ($("#searchResultsContainer").html().trim()!="No data found!") {
		if ($("#searchValue").val()){
			field = $("input[name=field]:checked").val()
			value = $("#searchValue").val()
			//Make sure the query was not already searched.
			if (current_query){
				for (var i in current_query) {
					if (current_query[i]["field"] == field && current_query[i]["value"] == value) {
						showRibbon("Already queried!", "#FFC87C","#popupContainer_inner", true);
						return false //Don't continue if query is already present.
					}
				}
			}

			showRibbon("Searching!", "#99FF5E","#popupContainer_inner", true);

			current_query.push({
				"field":field,
				"value":value,
			});
			sendSearchQuery(current_query);
		} else {
			showRibbon("Nothing entered!", "#FF6870", $("#popupContainer_inner"), true);
		}
	} else {
		showRibbon("No data to filter!", "#FF6870", $("#popupContainer_inner"), true);
	}
});

$(document).on("click", "#search_backButton", function() {
	if (current_query.length){
		current_query.pop();
		if (current_query.length){
			sendSearchQuery(current_query);
		} else {
			$("#searchResultsOuterContainer").html("<div id=\"searchResultsContainer\"> Enter filters to search. </div>");
		}
	} else {
		showRibbon("Already at start!", "#FF6870", $("#popupContainer_inner"), true);
	}
});

$(document).on("click", "#search_clearButton", function() {
	current_query = [];
	showRibbon("Filters emptied!", "#99FF5E","#popupContainer_inner", true);
	$("#searchResultsContainer").html("Enter filters to search.");
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
	if (activatorID == "leftMenu_addNew_copy") {
		activatorID = "leftMenu_addNew"
		event.stopPropagation();
	}

	$("#popupContainer").attr("for", activatorID);
	switch (activatorID) {
		case "leftMenu_addNew":
			var index = "";
			//Get the data to copy if applicable.
			if ($(this).attr("class").indexOf("duplicateSpecificDataButton") != -1) {
				index = $(this).siblings(".dataIndex").html().trim();
			}

			//Send the request to the server
			$.get("/data_form/"+index, function(response) {
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
			current_query = Array();
			$.get("/search/", function(response) {
				$("#popupContainer_inner").html(response);
				$("#tabs").tabs({active: 1});
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
		case "leftMenu_downloadCSV"://###Replace with cute form?
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

	//Remove any extra additions the popup may have populated.
	$(".CG_saveButton").remove()

	//Reload the screen if requested.
	if ($(".reloadActivator").length) {
		window.location.reload(true);
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

//######################################################################
});
