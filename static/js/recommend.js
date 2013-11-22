
$(document).ready(function() {
//############ Variable Setup: #########################################

$(document).on("mouseover", ".recommended_reactant", function() {
 if (CGEntries == undefined) {
  var specificDiv = $(this);
  $(specificDiv).attr("title", "Loading!")
  $.get("/send_CG_names/", function(response) {
   CGEntries = response;
   var compound = CGEntries[$(specificDiv).html().trim()] || "Compound not in guide!"
   $(".ui-tooltip").html(compound);
   });
 } else {
  var compound = CGEntries[$(this).html().trim()] || "Compound not in guide!"
  $(this).attr("title", compound)
 }
});

//######################################################################
});

