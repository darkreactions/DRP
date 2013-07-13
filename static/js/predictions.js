$(document).ready(function() {
//######################################################################
//Variable Setup
var dragging = false;
var dragStartX = 0;
var dragStartY = 0;
var widthSVG = 0;
var heightSVG = 0;
var dragResistance = 0.50;


//Load the SVG
var svg = $('#svgFile').svg(); 
svg.load(USER_SVG, {});
 
//Remove the title attribute on the first mouseover.
$(document).one("mouseover", "svg", function() {
	//Erase the SVG title.
	document.getElementById("graph1").getElementsByTagName("title")[0]
		.childNodes[0].nodeValue="";
	//Set the SVG size.
	widthSVG = document.getElementsByTagName("svg")[0].getAttribute("width");
	heightSVG= document.getElementsByTagName("svg")[0].getAttribute("height");
});

$(document).on("mousedown", "#svgFile", function(e) {
	dragging=true;
	dragStartX = e.pageX;
	dragStartY = e.pageY;
});
$(document).on("mouseup", "#svgFile", function(e) {
	dragging=false;
});

$("#svgFile").mousemove(function(e) {
	if (dragging) {
		dragEndX = e.pageX;
		dragEndY = e.pageY;
		
		
		var currentPos = document.getElementsByTagName("svg")[0].getAttribute("viewBox").split(" ");
		var newPos = []
		newPos.push(parseInt(currentPos[0]-(dragEndX-dragStartX)*dragResistance));
		newPos.push(parseInt(currentPos[1]-(dragEndY-dragStartY)*dragResistance));
		newPos.push(parseInt(currentPos[2]-(dragEndX-dragStartX)*dragResistance)+parseInt(widthSVG));
		newPos.push(parseInt(currentPos[3]-(dragEndY-dragStartY)*dragResistance)+parseInt(heightSVG));
		
		finalPos = newPos.join(" ");
		document.getElementsByTagName("svg")[0].setAttribute("viewBox", finalPos);
		
		dragStartX = e.pageX;
		dragStartY = e.pageY;
	}
});

//Node helper functions.
function nodeTooltip(node) {
	//Strip the title from the node.
	var elemID = $(node).attr("id");
	var elem = document.getElementById(elemID);
	$("#nodeTooltipContainer").html(
		elem.getElementsByTagName("title")[0].childNodes[0].nodeValue
		);
	elem.getElementsByTagName("title")[0].childNodes[0].nodeValue = ""
	
	//Get positioning variables.
	var nodePos = $(node).offset()
	var nodeCenterX = parseInt(nodePos.left) + parseInt($(node).children("ellipse").attr("rx"));
	var tooltipWidth = $("#nodeTooltipContainer").width();
	
	$("#nodeTooltipContainer").css({
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
	$("#nodeTooltipContainer").animate({
		opacity: 0,
	}, 200);
	
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
