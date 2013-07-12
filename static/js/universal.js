$(document).ready(function() {
//######################################################################
//############   Sorting   #############################################
//Sort from greatest to least.
window.sortNumbers = function(smallNum, bigNum) {
	return smallNum - bigNum;
}

window.sortNumbersReverse = function(smallNum, bigNum) {
	return bigNum - smallNum;
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

//Disable the tooltips when an input is selected.
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
	//alert($(".ribbonMessage").length);
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
	$.get("/user_logout/", function() {});
	
	//Reload the screen to verify log-out.
	refreshScreen()
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

$(document).on("click", "#search_filterButton", function() {
	if ($("#searchValue").val()){
		current_query.push({
			"field":$("input[name=field]:checked").val(),
			"value":$("#searchValue").val(),
		});
		alert(current_query);
		$.ajax({
			url:"/search/", 
			method:"post",
			data: {"current_query":JSON.stringify(current_query)}, 
			traditional: true,
			success: function(response) {
				$("#searchResultsOuterContainer").html(response)
			}
		});
	} else {
		showRibbon("Nothing entered!", "#FF6870", $("#popupContainer_inner"), true);
	}
	//WORKS FOR ONE QUERY:
	//if ($("#searchValue").val()){
		//current_query = {
			//"field":$("input[name=field]:checked").val(),
			//"value":$("#searchValue").val(),
		//};
		//$.post("/search/", current_query, function(response) {
			//$("#searchResultsOuterContainer").html(response)
		//});
	//} else {
		//showRibbon("Nothing entered!", "#FF6870", $("#popupContainer_inner"), true);
	//}
});

$(document).on("click", "#search_clearButton", function() {
	current_query = [];
});


//############ Popup Management: #######################################
window.createPopupConfirmation = function(message) {
	$("body").append("<div class=popupConfirmation title=Confirm:>"+
		message+"</div>");
}

//Activate specified popup:
$(document).on("click", ".popupActivator", function() {
	//Display a loading message if the request takes a visible amount of time.
	$("#popupContainer_inner").html("<center>Loading. Please wait.</center>");
	
	//Load the activator CSS.
	activatorID = $(this).attr("id");
	$("#popupContainer").attr("for", activatorID);
	switch ($(this).attr("id")) {
		case "leftMenu_addNew":
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
		case "leftMenu_uploadCSV":
			$.get("/upload_CSV/", function(response) {
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
alert("UNIVERSAL!");//###
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
