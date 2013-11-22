$(document).ready(function() {
//######################################################################
//Variable Setup
var dragging = false;
var dragStartX = 0;
var dragStartY = 0;
//Set the SVG size.
var widthSVG = parseInt(document.getElementsByTagName("svg")[1].getAttribute("width"));
var heightSVG= parseInt(document.getElementsByTagName("svg")[1].getAttribute("height"));
var dragResistance = 0.5;

var explorePath = [];
var nodeHoverGrowthProportion = 2;
var nodeAltered = false;


$(document).on("mousedown", "svg", function(e) {
 $(this).css("cursor","move");
 e.originalEvent.preventDefault();
 dragging=true;
 dragStartX = parseInt(e.pageX);
 dragStartY = parseInt(e.pageY);
});
$(document).on("mouseup", "svg", function(e) {
 $(this).css("cursor","pointer");
 dragging=false;
});

$(document).on("mouseleave", "#dataContainer", function(e) {
 $("svg").mouseup();
});

$(document).on("mousemove", "svg", function(e) {
 if (dragging) {
  dragEndX = parseInt(e.pageX);
  dragEndY = parseInt(e.pageY);
  
  var currentPos = $(this).get(0).getAttribute("viewBox").split(" ");
  
  //Make sure the SVG boundaries are not exceeded.
  var leftX = parseInt(currentPos[0])-parseInt(dragEndX-dragStartX)*parseFloat(dragResistance);
  if (leftX < 0) {leftX = 0}
  var topY = parseInt(currentPos[1])-parseInt(dragEndY-dragStartY)*parseFloat(dragResistance);
  if (topY < 0) {topY = 0}
   
  //Replace the current viewbox.
  var finalPos = [leftX, topY, widthSVG, heightSVG].join(" ");
  $(this).get(0).setAttribute("viewBox", finalPos);
  
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
  node.get(0).getElementsByTagName("title")[0].childNodes[0].nodeValue.substr(1).replace(/0/g," + ").replace(/_/g," ")
  );
 node.get(0).getElementsByTagName("title")[0].childNodes[0].nodeValue = ""
 
 //Get positioning variables.
 var nodePos = $(node).offset()
 var nodeCenterX = parseInt(nodePos.left) + parseInt($(node).children("ellipse").attr("rx"));
 var tooltipWidth = $("#nodeTooltipContainer").width();
 
 //Don't let the tooltips exceed the boundaries.
 var leftPos = nodeCenterX-tooltipWidth/2;
 if (leftPos<0) {
  leftPos = 0;
 }
 var upperPos = nodePos.top+30;
 if (upperPos + 60 + $("#nodeTooltipContainer").height() > window.innerHeight) {
  upperPos = nodePos.top-70 - $("#nodeTooltipContainer").height();
 }
 
 $("#nodeTooltipContainer").css({
  "display": "inline-block",
  "left": leftPos,
  "top": upperPos,
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
 if (!dragging) {
  nodeAltered = true;
  //Apply the nodeTooltip.
  nodeTooltip($(this));
  //Expand the node an itty-bit.
  var activeNode = $(this).get(0).getElementsByTagName("ellipse")[0];
  var radiusX = activeNode.getAttribute("rx");
  var newRadius = parseFloat(radiusX)*nodeHoverGrowthProportion;
  activeNode.setAttribute("rx", newRadius);
  activeNode.setAttribute("ry", newRadius);
 }
});

$(document).on("mouseover", ".edge", function() {
 //Erase the edge title on hover if a title is present..
 if ($(this).get(0).getElementsByTagName("title")[0].childNodes[0].nodeValue) {
  $(this).get(0).getElementsByTagName("title")[0].childNodes[0].nodeValue = "";
 }
});

$(document).on("mouseout", ".node", function() {
 $("#nodeTooltipContainer").fadeOut(200);
 
 //Reinsert the title node for future use.
 $(this).get(0).getElementsByTagName("title")[0].childNodes[0].nodeValue = $("#nodeTooltipContainer").html();

 if (nodeAltered) {
  //Shrink the node an itty-bit.
  var activeNode = $(this).get(0).getElementsByTagName("ellipse")[0];
  var radiusX = activeNode.getAttribute("rx");
  var newRadius = parseFloat(radiusX)/nodeHoverGrowthProportion;
  activeNode.setAttribute("rx", newRadius);
  activeNode.setAttribute("ry", newRadius);
  nodeAltered = false;
 }
});

//Node Clicks.
$(document).on("click", ".node", function() {
 //pastStep = document.getElementById("graph").getAttribute("class").split(" ")[1]
 //if (pastStep != "details") { //Details should be the last step.
  
  
  ////###The "Back" feature is contained below.
  //var source = $("#nodeTooltipContainer").html().toLowerCase().replace(/\+/g,"").replace(/ /g, "_").replace(/__/g,"_");
  var source = $("#nodeTooltipContainer").html().replace(/\+/g,"0").replace(/ /g, ""); //Make link-able.
  
  //###Jump over to Ironwood in a new window.
  var url = "http://ironwood.fig.haverford.edu/prefix/prefix_k"+source+".svg" 
  window.open(url, "_blank");
 
  ////Remember the SVG source/step.
  //explorePath.push([pastStep, source])
  //$.post("/gather_SVG/", {
    //"step":pastStep,
    //"source":source
   //}, 
   //function(response){
    //$("#dataContainer").html(response);
    
    //if ($(".fatalError").length>0){
     //$("#nodeTooltipContainer").remove()
     //return false;
    //} 
    
    ////Set the SVG title to the new "graph step."
    //createSVGTitle(document.getElementById("graph").getAttribute("class").split(" ")[1]);
    //createSVGBack();
  //})
  //return false; //Don't continue or the graph will be clicked.
 //}
});

//#######################   SVG Buttons    #############################
$(document).on("click", ".svgBackButton", function() {
 //Take the last graph in the explorePath.
 explorePath.pop();
 if (explorePath.length>0) {
  var lastGraph = explorePath[explorePath.length-1]
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
   createSVGTitle($("#svgControl").attr("currentStep"));
   if (lastGraph[0]!="start") {
    createSVGBack();
   }
 })
 return false; //Don't continue or the graph will be clicked.
});

//#######################   On Load    #################################
//Set up the SVG variables once everything has loaded.
$("svg").css("cursor","pointer");


//Finally, erase the main SVG title.
document.getElementById("graph1").getElementsByTagName("title")[0].childNodes[0].nodeValue="";

//Create the SVG Title.
svgControl = "<div id=\"svgControl\" currentStep=\"compounds\"></div>"
$("body").append(svgControl);
createSVGTitle($("#svgControl").attr("currentStep"));


//######################################################################
});
