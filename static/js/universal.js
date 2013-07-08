$(document).ready(function() {
//######################################################################
//############   Window Management   ###################################
//Set the current_page cookie.
function setPageCookie(page) {
	alert("3");
	page = page !== undefined ? page : $("#pagesCurrent").html().trim();
	$.cookie("current_page", page,
		{expires: 7}); //Set cookie to expire after one week.
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
//############ User Authentication: ####################################
$("#userLogOut").click( function() {
	//Send the log-out signal.
	$.get("/user_logout/", function() {});
	
	//Reload the screen to verify log-out.
	refreshScreen()
});

//############ Popup Management: #######################################
function createPopupConfirmation(message) {
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

//######################################################################
});
