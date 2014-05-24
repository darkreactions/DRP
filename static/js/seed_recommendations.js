$(document).ready(function(){
//######################################################################

var timer = applyLoadingWheel("#seedRecOven");

//Every thirty seconds, check the seed-rec oven.
setInterval(function() {
  //Clear the last "timer" so that two aren't set simultaneously.
  clearInterval(timer);

  $.get("/check_seed_oven/", function(response){
    $("#seedRecOven").replaceWith(response);
    timer = applyLoadingWheel("#seedRecOven");
  });
}, 6*1000);


//######################################################################
});
