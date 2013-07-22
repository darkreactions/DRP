$(document).ready(function() {
//######################################################################
//Variable Setup
var dragging = false;
var dragStartX = 0;
var dragStartY = 0;
var widthSVG = 0;
var heightSVG = 0;
var dragResistance = 0.5;

var explorePath = []

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

$(document).on("mousemove", "#svgFile", function(e) {
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

function createSVGTitle(text) {
	$(".svgTitle").remove();
	var titleContainer = "<div class=\"svgTitleContainer svgTextUI\">" + text[0].toUpperCase() + text.substr(1) + "</div>" 
	$("#dataContainer").append(titleContainer);
}

function createSVGBack() {
	$(".svgBack").remove();
	var titleContainer = "<div class=\"svgBackButton genericButton svgTextUI\">Back</div>" 
	$("#dataContainer").append(titleContainer);
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
	pastStep = document.getElementById("graph").getAttribute("class").split(" ")[1]
	if (pastStep != "details") { //Details should be the last step.
		
		var source = $("#nodeTooltipContainer").html().toLowerCase().replace(/\+/g,"").replace(/ /g, "_").replace(/__/g,"_");
	
		//Remember the SVG source/step.
		explorePath.push([pastStep, source])
		$.post("/gather_SVG/", {
				"step":pastStep,
				"source":source
			}, 
			function(response){
				$("#dataContainer").html(response);
				
				if ($(".fatalError").length>0){
					$("#nodeTooltipContainer").remove()
					return false;
				} 
				
				//Set the SVG title to the new "graph step."
				createSVGTitle(document.getElementById("graph").getAttribute("class").split(" ")[1]);
				createSVGBack();
		})
		return false; //Don't continue or the graph will be clicked.
	}
});

//#######################   SVG Buttons    #############################
$(document).on("click", ".svgBackButton", function() {
	//Take the last graph in the explorePath.
	explorePath.pop();
	if (explorePath.length>0) {
		var lastGraph = explorePath.pop()
	} else {
		var lastGraph = ["start", null];
	}
	
	$.post("/gather_SVG/", {
			"step":lastGraph[0],
			"source":lastGraph[1]
		}, 
		function(response){
			$("#dataContainer").html(response);
			
			if ($(".fatalError").length>0){
				$("#nodeTooltipContainer").remove()
				return false;
			} 
			
			//Set the SVG title to the new "graph step."
			createSVGTitle(document.getElementById("graph").getAttribute("class").split(" ")[1]);
			if (lastGraph[0]!="start") {
				createSVGBack();
			}
	})
	return false; //Don't continue or the graph will be clicked.
});

//#######################   On Load    #################################
createSVGTitle(document.getElementById("graph").getAttribute("class").split(" ")[1]);

//######################################################################
});
