$(document).ready(function() {
//######################################################################
//Variable Setup
var dragging = false;
var dragStartX = 0;
var dragStartY = 0;
var widthSVG = 0;
var heightSVG = 0;
var dragResistance = 0.5;

$(document).on("click", "#svgReload", function() {
	showRibbon("Trying!", "#99FF5E", "#dataContainer");
	loadSVG();
});


//Remove the title attribute on the first mouseover.
$(document).one("mouseover", "#svgFile", function() {
	//Set the SVG size.
	widthSVG = parseInt(document.getElementsByTagName("svg")[0].getAttribute("width"));
	heightSVG= parseInt(document.getElementsByTagName("svg")[0].getAttribute("height"));
	
	//Finally, erase the SVG title.### Needed?
	//document.getElementById("graph1").getElementsByTagName("title")[0]
	//	.childNodes[0].nodeValue="";
});

$(document).on("mousedown", "#svgFile", function(e) {
	$("#svgFile").css("cursor","move");
	e.originalEvent.preventDefault();
	dragging=true;
	dragStartX = parseInt(e.pageX);
	dragStartY = parseInt(e.pageY);
});
$(document).on("mouseup", "#svgFile", function(e) {
	$("#svgFile").css("cursor","pointer");
	dragging=false;
});

$(document).on("mouseout", "#svgFile", function(e) {
	$("#svgFile").css("cursor","pointer");
	dragging=false;
});

//////$(document).on("mousewheel", "#svgFILE", function(e) {###
	//////alert("cat");
//////});

$("#svgFile").mousemove(function(e) {
	if (dragging) {
		dragEndX = parseInt(e.pageX);
		dragEndY = parseInt(e.pageY);
		
		
		var currentPos = document.getElementsByTagName("svg")[0].getAttribute("viewBox").split(" ");
		
		//Make sure the SVG boundaries are not exceeded.
		var leftX = parseInt(currentPos[0])-parseInt(dragEndX-dragStartX)*parseFloat(dragResistance);
		if (leftX < 0) {leftX = 0}
		var topY = parseInt(currentPos[1])-parseInt(dragEndY-dragStartY)*parseFloat(dragResistance);
		if (topY < 0) {topY = 0}
		 
		//Replace the current viewbox.
		var finalPos = [leftX, topY, widthSVG, heightSVG].join(" ");
		document.getElementsByTagName("svg")[0].setAttribute("viewBox", finalPos);
		
		dragStartX = dragEndX;
		dragStartY = dragEndY;
	}
});

//Node helper functions.
function nodeTooltip(node) {
	$("#nodeTooltipContainer").remove();
	
	//Strip the title from the node.
	$("body").append("<div id=\"nodeTooltipContainer\"></div>");
	$("#nodeTooltipContainer").html(
		node.get(0).getElementsByTagName("title")[0].childNodes[0].nodeValue
		);
	node.get(0).getElementsByTagName("title")[0].childNodes[0].nodeValue = ""
	
	//Get positioning variables.
	var nodePos = $(node).offset()
	var nodeCenterX = parseInt(nodePos.left) + parseInt($(node).children("ellipse").attr("rx"));
	var tooltipWidth = $("#nodeTooltipContainer").width();
	
	$("#nodeTooltipContainer").css({
		"display": "inline-block",
		"left": nodeCenterX-tooltipWidth/2,
		"top": nodePos.top+50,
	});
	
	$("#nodeTooltipContainer").animate({
		opacity: 1,
	}, 500);
}

function nodePopup(node, startX, startY) {
	
}

//Node Mouseovers
$(document).on("mouseover", ".node", function() {
	nodeTooltip($(this));
});
$(document).on("mouseout", ".node", function() {
	$("#nodeTooltipContainer").fadeOut(200);
	
	//Reinsert the title node for future use.
	$(this).get(0).getElementsByTagName("title")[0].childNodes[0].nodeValue = $("#nodeTooltipContainer").html();
});

//Node Clicks.
$(document).on("click", ".node", function() {
	alert("Do something! : )");
	return false; //Don't continue or the graph will be clicked.
});
 
//######################################################################
});
