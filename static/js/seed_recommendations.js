$(document).ready(function(){
//######################################################################

//Every thirty seconds, check the seed-rec oven.
setInterval(function() {
  $.get("/check_seed_oven/", function(response){
    $("#seedRecOven").replaceWith(response);
    applyLoadingWheel("#seedRecOven");
  });
}, 6*1000);

applyLoadingWheel("#seedRecOven");

//######################################################################
});
