$(document).ready(function() {
//######################################################################
//Variable Setup
var dragging = false;
var dragStartX = 0;
var dragStartY = 0;
var widthSVG = 0;
var heightSVG = 0;
var dragResistance = 0.5;

//Load the SVG
window.loadSVG = function() {
	svg = $('#svgFile').svg(); 
	svg.load(USER_SVG, {});
}
loadSVG();

$(document).on("click", "#svgReload", function() {
	showRibbon("Trying!", "#99FF5E", "#dataContainer");
	loadSVG();
});


//Remove the title attribute on the first mouseover.
$(document).one("mouseover", "svg", function() {
	//Erase the SVG title.
	document.getElementById("graph1").getElementsByTagName("title")[0]
		.childNodes[0].nodeValue="";
	//Set the SVG size.
	widthSVG = parseInt(document.getElementsByTagName("svg")[0].getAttribute("width"));
	heightSVG= parseInt(document.getElementsByTagName("svg")[0].getAttribute("height"));
	$("#svgFile").css("cursor","pointer");
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

$(document).on("mousewheel", "#svgFILE", function(e) {
	alert("cat");
});

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
	var elemID = $(node).attr("id");
	var elem = document.getElementById(elemID);
	$("body").append("<div id=\"nodeTooltipContainer\"></div>");
	$("#nodeTooltipContainer").html(
		elem.getElementsByTagName("title")[0].childNodes[0].nodeValue
		);
	elem.getElementsByTagName("title")[0].childNodes[0].nodeValue = ""
	
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
	//Change the node color.
	$(this).children("ellipse").attr("stroke-width" ,"2");
	$(this).children("ellipse").attr("stroke", "yellow");
	nodeTooltip($(this));
});
$(document).on("mouseout", ".node", function() {
	$(this).children("ellipse").attr("stroke-width" ,"1");
	$(this).children("ellipse").attr("stroke", "black");
	//$("#nodeTooltipContainer").animate({
		//opacity: 0,
		//zindex:0,
	//}, 200);
	$("#nodeTooltipContainer").fadeOut(200);
	
	//Reinsert the title node for future use.
	var elem = document.getElementById($(this).attr("id"));
	elem.getElementsByTagName("title")[0].childNodes[0].nodeValue = $("#nodeTooltipContainer").html();
});

//Node Clicks.
$(document).on("click", ".node", function() {
	alert("Do something! : )");
	return false; //Don't continue or the graph will be clicked.
});
 
//######################################################################
});
