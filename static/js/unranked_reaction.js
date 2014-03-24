$(document).on("ready", function() {

$(".sortable").sortable();
$(".sortable").disableSelection();

$("#rankSubmit").click(function() {
	if (window.confirm("Is this your final answer?")==false) {
		return false;
	} else if ($(".unrankedReaction").length==0){
   		showRibbon("No reactions to send!", badColor, "#mainPanel");
		return false;
	} else {
   		showRibbon("Loading!", neutralColor, "#mainPanel");
	}

	//Get the new order.
	var newOrder = Array();
	$.each( $(".unrankedReaction.sorted "), function(){
		var reactionList = Array();
		$.each( $(this).find(".dataField"), function(){
			reactionList.push($(this).html().trim());
		});
		newOrder.push(reactionList);
 	})

	$.post("/send_and_receive_rank/",
		{
			"pid": $(".sortable").attr("pid"), 
			"newOrder":JSON.stringify(newOrder)
		}, function(response) {
		 	$("#rankingContainer").html(response)	
			$(".sortable").sortable();
			$(".sortable").disableSelection();
		});

});

})

