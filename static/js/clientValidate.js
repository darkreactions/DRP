//Load necessary JSON files.
var dataRanges;
$.ajax({
	"async": false,
	"dataType": "json",
	"url":STATIC_URL+"/js/dataRanges.json",
	"success": function(json) {
		dataRanges = json;
	}
});

var editChoices;
$.ajax({
	"async": false,
	"dataType": "json",
	"url":STATIC_URL+"/js/editChoices.json",
	"success": function(json) {
		editChoices = json;
	}
});
	
//Fields that can be edited with a range alone: ###Copied in validation.py
var rangeFields = ["quantity", "temp", "time", "pH", "outcome", "purity"];
//Fields that must be between a specific character count/limit:
var limitFields = ["ref", "notes"];

function quickValidate(field, value) {//###
	//Get the range if it is applicable.
	if (dataRanges[field]) {
		var range = dataRanges[field];
	}
	
	if (rangeFields.indexOf(field) >= 0){ 
		return (range[0] <= parseFloat(value) && parseFloat(value) <= range[1])
	} else if (limitFields.indexOf(field) >= 0) {
		return (range[0] <= value.length && value.length <= range[1])
	} else {
		return true; //If range tests passed, label as valid data.
	}
}

function fullValidate(field, value) {
	//Name Validation
	
	return quickValidate(field, value); //Finish with range validations.
}
