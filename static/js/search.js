
$(document).on("keyup", "#searchValue", function(e) {
  if (e.keyCode==13) { // The `enter` key
    $(".search_filterButton").trigger("click");
  }
})

$(document).on("click", ".search_resetFiltersButton", function() {
  location.search = "";
});

$(document).on("click", ".search_filterButton", function() {

  var $form = $(this).closest("form");
  var field = $form.find("select[name='field']").val();
  var subfield = $form.find("select[name='subfield']:visible").val();
  var match = $form.find("select[name='match']").val();
  var value = $form.find("#searchValue").val();

  if (value!=="" && value!==undefined) {

    var key = field;
    if (subfield!==undefined) key += "."+subfield;
    if (match!==undefined) key += "."+match;

    // Set the new querystring and reload the page.
    var newQueryParam = "&" + window.encodeURI(key+"="+value)
    location.search += newQueryParam;

  } else {
    var failureMessage = "No query entered!"
    showRibbon(failureMessage, badColor, "#mainPanel");
  }

});



// Set the autocomplete box on reactants if needed.
$(document).on("change", ".dropDownMenu[name=field]", function() {
  if ($(this).val()=="reactant"){
    $(".autocomplete_reactant").autocomplete("enable")
    $("#searchValue").addClass("autocomplete_reactant")
  } else {
    $(".autocomplete_reactant").autocomplete("disable")
  }
});

// Disable the "subfield" if it is unnecessary.
$(document).on("change", "#searchForm .dropDownMenu[name=field]", function(){
  var val = $(this).find("option:selected").val();

  if (val==="reactant") {
    $("#searchForm select[name=subfield]").show();
  } else {
    $("#searchForm select[name=subfield]").hide();
  }
});

