$(document).ready(function() {

//Set the current_page cookie.
function setPageCookie(page, url) {
 page = page !== undefined ? page : $("#pagesCurrent").html().trim();
 $.cookie(url+"current_page", page,
  {expires: 7}); //Set cookie to expire after one week.
}


function changePageTo(page) {
 showRibbon("Loading!", "#FFC87C", "#mainPanel", false);

 //Send the currentQuery if it exists.
 JSONQuery = JSON.stringify({"currentQuery":currentQuery, "page":page});

 var attribute = $(".pagesInfo").attr("model");
 var url = (attribute==="database") ? "/data_transmit/" : "/recommend_transmit/";

 $.post(url, {"body":JSONQuery}, function(response) {
  $("#mainPanel").html(response);
  return;
  if ($(".dataGroup").length) {
   setPageCookie(page, url)
   restyleData();
  $(".ribbonMessage").remove();
  } else {
   showRibbon("Page does not exist", badColor, "body");
  }
 });
}

$(document).on("click", "#pagesInputButton", function() {
 var pageNum = $("#pagesInputText").val().trim();

 if (pageNum == $("#pageLinkCurrent").html().trim()) {
  showRibbon("Already here!", "#99FF5E", "#dataContainer");
  return false //Don't do anything else if page already is active.
 } else if ($.isNumeric(pageNum)
  && parseInt(pageNum) <= parseInt($("#pagesTotal").html().trim())
  && parseInt(pageNum) > 0) {
  //Create a CSRF Token for Django authorization.
  var csrftoken = $.cookie("csrftoken");
  changePageTo(pageNum)
  }
});

$(document).on("click", ".pageLink", function() {
 //If a valid page link has been clicked.
 var pageNum = $(this).html().trim();
 //If the link is a number (not an ellipsis) and not the current page:
 if ($.isNumeric(pageNum) && pageNum != $("#pageLinkCurrent").html().trim()) {
  //Request the information about a page.
  changePageTo(pageNum)
 }
});


//######################################################################
});
