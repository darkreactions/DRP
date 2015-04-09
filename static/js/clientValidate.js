//Load necessary JSON files.
var dataRanges;
$.ajax({
 "async": true,
 "dataType": "json",
 "url":STATIC_URL+"/js/dataRanges.json",
 "success": function(json) {
  dataRanges = json;
 }
});

var editChoices;
$.ajax({
 "async": true,
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
var notRequired = ["CAS_ID", "pH"]

function quickValidate(field, value, required) {//###
 required = required !== undefined ? required : true;
 if ((!required && value=="") || (notRequired.indexOf(field)>=0 && value=="")) {
  return true;
 }

 //Get the range if it is applicable.
 if (dataRanges[field]) {
  var range = dataRanges[field];
 }

 if (rangeFields.indexOf(field) >= 0){
  return (range[0] <= parseFloat(value) && parseFloat(value) <= range[1])
 } else if (limitFields.indexOf(field) >= 0) {
  return (range[0] <= value.length && value.length <= range[1])
 } else if (field=="reactant"){
  //Load the CG entries if not loaded yet.
  if (CGEntries == undefined) {
   $.get("/send_CG_names/", function(response) {
    CGEntries = response;
    });
  }

  if (CGEntries[value]!==undefined) {
   return true;
  }
  return false;
 } else {
  return true; //If range tests passed, label as valid data.
 }
}
