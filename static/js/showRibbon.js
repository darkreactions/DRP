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
   $(location).children(".ribbonMessage").fadeOut(1000, function()
   {
    $(".ribbonMessage").remove();
   });

  },500+15*message.length);
 }
}
