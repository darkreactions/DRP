//All range attributes can be modified CLIENT-SIDE by changing these ranges.
//Note: ranges were arbitrarily chosen for unit protection and should
//	be adjusted. 	
// (To change the ranges SERVER-SIDE, modify data_ranges.py.)
var quantity_range = [0,15]; //grams
var temp_range = [0,400]; //Celsius
var outcome_range = [0,4]; //0=bad, 4=good
var purity_range = [0,2]; //0=bad, 2=good
var time_range = [0,350]; //hours
var pH_range = [-1,16]; //[unit-less]
var ref_range = [0,8]; //Number of characters in the string.

var massFields = ["quantity_1","quantity_2", "quantity_3", "quantity_4", "quantity_5"];
var dropsPermL = 20;
var mLConversions = { //Must be strings to avoid simplification by Javascript.
	"EtOH":"0.789", //Units must be g/mL
	"H2O":"1.000", 
}


//Click-to-edit Data
var editByMenu = { //###FIX ME
	"type_outcome":[0,1,2,3,4],
	"type_purity": [0,1,2],
	"type_unit":["g", "mL", "d"],
	"type_slow_cool":["Yes", "No", "?"], //###
	"type_leak":["Yes", "No", "?"]
};

$(document).on("click", ".sub_data", function() {
	if ($(this).children().is(":focus")) {
		//If the element is already in focus, don't do anything.
	} else {
		var originalText = $(this).text().trim();
		var originalWidth = $(this).css("width");
		var data_type = $(this).attr("class").split(' ')[1];
		
		// EDIT BY DROP-DOWN MENU.
		if (editByMenu.hasOwnProperty(data_type)) {
			var new_innards = ("<select id=\"editInPlaceMenu\" "
				+"originalText=\"" + originalText + "\">"
			);
			//Construct the edit menu with the current value selected.
			options = editByMenu[data_type]
			for (var i in options) {
				if (options[i] == originalText) { //Auto-select the current value.
					new_innards += "<option value=" + options[i] + " selected>" + options[i] + "</option>";
				} else {
					new_innards += "<option value=" + options[i] + ">" + options[i] + "</option>";
				}
			}
			
			new_innards += "</select><div id=\"enterButton\">Enter</div>";
			$(this).html(new_innards);
			$(this).children().focus();
			
		} else {
			$(this).html(
				"<input id=\"editInPlace\" type=\"text\""
				+ " originalText=\"" + originalText + "\" "
				+ " value=\"" + originalText + "\" style=\"width:" 
				+ String(parseInt(originalWidth)) + ";\" />"
				+ "<div id=\"enterButton\"> Enter</div>"
			);
			$(this).children().focus();
		}
	}
});


//Replace #editInPlace elements with just the value when it loses focus.
$(document).on("blur", "#editInPlaceMenu", function() {
	submitEditToServer($(this));
});

//Replace #editInPlace elements with just the value when it loses focus.
$(document).on("blur", "#editInPlace", function() {
	submitEditToServer($(this));
});



//||||||||| Client-Side Validation: ||||||||||||||||||||||||||||||||||||
//Quickly test if data is the right type and in the correct range.
function quickValidation(field, dirty_datum) { //(Mainly ported from validation.py
	try {
		var dirty_datum;
		if (massFields.indexOf(field) != -1) {
			dirty_datum = Number(dirty_datum);
			return(
				(String(dirty_datum).length<10) && //#Prevent super long decimals.
				(typeof (dirty_datum) == typeof (1.0)) &&
				(quantity_range[0] < dirty_datum) && (dirty_datum  < quantity_range[1])
			);
		} else if (field=="temp") {
			dirty_datum = Number(dirty_datum);
			return (
				(temp_range[0] < dirty_datum) && (dirty_datum  < temp_range[1])
			);
		} else if (field=="pH") {
			dirty_datum = Number(dirty_datum);
			return (
				(pH_range[0] < dirty_datum) && (dirty_datum  < pH_range[1])
			);
		} else if (field=="outcome") {
			dirty_datum = Number(dirty_datum);
			return (
				(outcome_range[0] <= dirty_datum) && (dirty_datum  <= outcome_range[1])
			);
		} else if (field=="purity") {
			dirty_datum = Number(dirty_datum);
			return (
				(purity_range[0] <= dirty_datum) && (dirty_datum  <= purity_range[1])
			);
		} else if (field=="time") {
			dirty_datum = Number(dirty_datum);
			return (
				(time_range[0] < dirty_datum) && (dirty_datum  < time_range[1])
			);
		} else if (field=="slow_cool" || field=="leak") {
			dirty_datum = String(dirty_datum).toLowerCase()
			if (["yes","no","?"].indexOf(dirty_datum) != -1) { //UNKNOWN LABEL
				return (true);
			} else {
				return (false);
			} 
		} else if (field=="ref") {
			dirty_datum = String(dirty_datum)
			return ( 
				(ref_range[0] < dirty_datum.length) && 
				(ref_range[1] >= dirty_datum.length)
				);
		}
		return true; //Should not be reached if all fields are validated above.
	} catch (error) {
		return false; //If something goes wrong, assume false data.
	}	
}

$(".fieldWrapper input").blur(function() {
	//Prepare Variables
	var dataIsValid = true; //Boolean that changes in regard to data tests.
	var dirtyField = $(this).attr("id").substr(3);
	var dirtyData = $(this).val();
	
	//CLIENT-SIDE VALIDATION ###
	//Let the user continue if no data or valid data is entered.
	dataIsValid = quickValidation(dirtyField, dirtyData) || !dirtyData
	message = "Invalid Data"; //###Differentiate?
	
	if (dataIsValid) {
		//Remove "invalid" style if applicable.
		if ($(this).attr("class").split(" ").indexOf("invalid_data") != -1) {
			$(this).removeClass("invalid_data");
		}
	} else {
		//Display a message to the user to let them know their data is bad.
		displayErrorMessage(message);

		//Mark incorrect fields, but only mark them once.
		if ($(this).attr("class") != "invalid_data") {
			$(this).addClass("invalid_data");
		}
	}
});

//||||||||| User Authentication: ||||||||||||||||||||||||||||||||||||
$(".userAuthLink").click( function() {
	//Get the Registration form from the server if it is requested.
	if ($(this).attr("id") == "userRegister") {
		$.get("/user_registration/", function(formHTML) {
			$("#userForm").html(formHTML);
			$("#id_username").focus(); //###Should just be "first visible input"
		});
	} //Get the Log-In form if it is requested.
	else if ($(this).attr("id") == "userLogin") {
		$.get("/user_login/", function(formHTML) {
			$("#userForm").html(formHTML);
			$("#id_username").focus(); //###Should just be "first visible input"
		});
	}
	
	
	//Show and center the userAuthScreen.
	$("#userAuthScreen").css({
		"margin":"0px auto",
		"top":"0px",
		"left":"0px",
		"display":"block"});
	$("#popupGlobal").show();
});

$("#userLogOut").click( function() {
	//Send the log-out signal.
	$.get("/user_logout/", function(formHTML) {});
	
	//Reload the screen to verify log-out.
	//refreshScreen()
});

//||||||||| Upon Form Submissions: ||||||||||||||||||||||||||||||||||||


$(document).on("submit", "#dataEntryForm", function() { //###Not generalized to all popups.
	var form = $(this); //Keep a reference to the form inside the POST request.
	//$.post($(form).attr("action"), $(form).serialize(), function(response) {
		//Recreate the popup window with the server response.
		//alert(response);
	//});
	//alert("yes!");
	//return false; //Do not continue or else the form will post. 
});

//||||||||| Upon Page Load: ||||||||||||||||||||||||||||||||||||

//Load the error log if values were not valid server-side.
if (!$("#error_container_text").html().trim()=="") {
	$("#error_container").css({
		"margin":"0px auto",
		"position":"fixed",
		"top":"30%",
		"display":"block"});
		
	$("#error_container").show();
	$("#error_container").draggable();
}

$(document).on("click", "#popup_toggle", function() { //###Not generalized to all popups.
	$("#error_container_text").animate({
		height:"toggle",
		opacity:"toggle",
		},500, function(){});
});

$(document).on("click", "#popup_close", function() {
	var popupScreen = $(this).parent();
	if (popupScreen.parent().attr("id") == "popupGlobal") {
		popupScreen.parent().animate({
			height:"toggle",
			opacity:"toggle",
			},500, function(){});
	} else {
		popupScreen.animate({
			height:"toggle",
			opacity:"toggle",
			},500, function(){});
	}
});

refreshDataContainer();

});//Close the ready() call at document head.
