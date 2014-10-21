//var nodePositionsPath = "/home/ntien/workspace/dark-reaction-site/DRP/vis/nodePositions.json" 
var width = parseInt(d3.select("#graph").style("width"), 10);
var height = parseInt(d3.select("#graph").style("height"), 10); 
d3.json("/get_graph/", function(graph) {
    
    var nodeTooltips = [["Purity", "purity"],
			["Outcome", "outcome"],
			["Reference Number", "ref"], 
			["Inorg1", "inorg1"],
			]; 
    
    var nodes = graph.nodes;   
    var preLoad = graph.skipTicks === "True";  
    var links = graph.links; 
    var clusters = graph.clusters;
    //declaring a variable that will serve to size the text boxes that appear upon hovering over the individual nodes
    var textPlacement = nodeTooltips.length*40;
    var boxLength = 130;
    

  var translateValue1 = 500;
  var translateValue2 = 300;
  var scaleValue =.2; 

  var zoom = d3.behavior.zoom()
    .scaleExtent([1/10, 2])
    .on("zoom", zoomed)
    .translate([translateValue1, translateValue2]).scale(scaleValue); 

  var drag = d3.behavior.drag()
      .origin(function(d) { return d; })
      .on("dragstart", dragstarted)
      .on("drag", dragged)
      .on("dragend", dragended);

//Needed for zooming and dragging (http://bl.ocks.org/mbostock/6123708).
function dragstarted(d) {
  d3.select(this).classed("dragging", true);
}

function dragged(d) {
  d3.select(this)
    .attr("cx", d.x = d3.event.x)
    .attr("cy", d.y = d3.event.y);
}

function dragended(d) {
  d3.select(this).classed("dragging", false);
}

    //Add a default weight of   1   for each connection.
    links = links.map(function(connection) {
      if (connection["source"]<0 || connection["target"]<0) {
        console.log(connection);
      }
      connection["weight"] = 1;
      return connection;
    })
  
var force = d3.layout.force()
  .nodes(nodes)
  .links(links)
  .charge(-500)
  .friction(0.6)
  .gravity(0.4)
  .linkDistance(100)
  .size([width, height]);
 
var svg = d3.select("#graph").append("svg")
    .attr("width", width)
    .attr("height", height)
    .call(zoom);

var centerX = width / 2 
var centerY = height / 2 

var container = svg.append("g")

container.append("rect")
	.attr("width", width)
	.attr("height", height)
	.attr("fill", "white") 

container.call(zoom).attr("transform","translate(" + translateValue1 +"," + translateValue2 +")scale(" + scaleValue +"," + scaleValue +")");


function zoomed() {
  container.attr("transform", 
  "translate(" + d3.event.translate +")scale("+ d3.event.scale + ")");
}

// Use a timeout to allow the rest of the page to load first.
setTimeout(function() {

  // Run the layout a fixed number of times.
  // The ideal number of times scales with graph complexity.
  // Of course, don't run too long—you'll hang the page!
var maxIterations = parseFloat(100); 
var loadingBar = $("#innerLoadingBar");
var loadingBarMaxLength = loadingBar.parent().width();


if (!preLoad) {
 console.log(preLoad) 
 console.log("made it to preload!")  
 force.on("tick", function() {
  nodes[0].x = width / 2;
  nodes[0].y = height / 2;
 }); 

  force.start();
  for (var i = 0; i < maxIterations; i++) {
	force.tick()
  }
  force.stop();  
  } else {
   force.start();  
   force.tick();
   force.stop();
  } 
    
    
 
  container.selectAll("line")
    .data(links)
  .enter().append("line")
    .attr("x1", function(d) { return d.source.x; })
    .attr("y1", function(d) { return d.source.y; })
    .attr("x2", function(d) { return d.target.x; })
    .attr("y2", function(d) { return d.target.y; })
    .style("stroke-width", 0.06)
    .attr("stroke", "gray")
/*
  var circleClusters = container.selectAll(".clusters")
      .data(clusters) 
      .enter()
      .append("circle") 
      .attr("cx", function(d) {return d.x;})
      .attr("cy", function(d) {return d.y;})
      .attr("r", "20")
      .attr("fill", function(d) {return d.color;})
      .attr("opacity", "0.5") 
*/ 
container.on("mouseover", function() {d3.selectAll(".tooltipContainer").remove(); console.log("made it!");}); 
 
    var nodeElements = container.selectAll(".node")
    .data(nodes.filter(function(d) { return d.outcome > 0;}))

    .enter().append("g")
    .attr("class", "node"); 

   
    nodeElements.append("circle").attr("class", ".circleClusters").attr("fill", function(d) {return (d.color!=undefined) ? d.color : "rgba(0,0,0,0)";}).attr("opacity", 0.4).attr("r", 20) 

    nodeElements.append("circle") 
    .attr("r", function(d) { var size = Math.abs(Math.log(d.pagerank))/3 + d.pagerank*450;
    if (size > 10){
      size = 10
      }
    else if (size < 0.5){
      size = 0.5;
      }
      return size
      })
    .style("fill", function(d) {if (d.outcome == 4){
                    return "#1a9641";
                    }
                  else if (d.outcome == 3){
                    return "#a6d96a";
                    }
                  else if (d.outcome == 2){
                    return "#fdae61"
                    }
                  else if (d.outcome == 1) {
                    return "#d7191c";
                    }
                  else {
                    return "purple";
                    }
                  })
    .on("click", function() {
	d3.event.stopPropagation()})

    .on("click", function() {
	d3.selectAll(".tooltipContainer").remove(); 
	var thisGroup = this.parentNode;
	this.parentNode.parentNode.appendChild(thisGroup); 
    	
	var currentCircle = d3.select(thisGroup);  
	
	var textbox = currentCircle.append("g")
	  .attr("class", "tooltipContainer") 
        
	textbox.append("rect")
	  .attr("class", "tooltipBackground")
	  .attr("width", 350)
	  .attr("height", boxLength * 1.3)
      .on("mouseover", function() {
	d3.event.stopPropagation()})

  
	var defs = container.append("defs"); 

	var filter = defs.append("filter") 
	  .attr("id", "drop-shadow")
	  .attr("height", "130%")
	  .attr("width", "130%"); 
	
	filter.append("feGaussianBlur")
	  .attr("in", "SourceAlpha")
	  .attr("stdDeviation", 1)
	  .attr("result", "blur"); 
	filter.append("feOffset")
	  .attr("in", "blur")
	  .attr("dx", 1)
	  .attr("dy", 1)
	  .attr("result", "offsetBlur");
	var feMerge = filter.append("feMerge");
	feMerge.append("feMergeNode")
	  .attr("in", "offsetBlur")
	feMerge.append("feMergeNode")
	  .attr("in", "SourceGraphic");

			
	var textElement = textbox.append("text")
	  .attr("class", "tooltip")
	
	var d = thisGroup.__data__
	for (var i=0; i < nodeTooltips.length; i++) { 

	  var fieldName = nodeTooltips[i][0] 
	  var fieldValue = d[nodeTooltips[i][1]]
	  var textField = fieldName + ": " + fieldValue 
 	  textElement.append("tspan")
		   .text(textField) 
		   .attr("x", "20px")
		   .attr("dy", "2em"); 
	  }; 
	var inorg2 = d["inorg2"];
	if(typeof inorg2 !== "undefined") 
	  {textElement.append("tspan")
		.text("Inorg2: " + inorg2)
		.attr("x", "20px")
		.attr("dy", "2em")
       boxLength = 160;
	} else {
    boxLength = 130; 
	}; 
    var seedRecButton = textbox.append("g")
		.attr('cursor', 'pointer')
		.on("mouseover", function() {
		seedRecButton.select("rect").style("filter","url(#drop-shadow)");}) 
		.on("mouseout", function() { 
		seedRecButton.select("rect").style("filter", "none");}) 

	seedRecButton.append("rect")
		.attr("class", "seedButtonRect")  
		.attr('width', 310)
		.attr('height', 25)  
		.attr('x', "20px")
		.attr('y', boxLength)
		.attr("rx", "3")
		.attr("ry", "3");
	//	.style("filter", "url(#drop-shadow)");
	
	seedRecButton.append("text")
		.text("Generate seed recommendations")
        .attr("class", "seedRecText")
		.attr('x', "25px")
		.attr('y', boxLength + 15) 
		.append("rect")
		.attr("fill","none"); 
 	
	var imageContainer = seedRecButton.append("g")
	var seedButton = imageContainer.append("svg:image")
	  .attr('x', textPlacement*1.7) 
 	  .attr('y', boxLength + 5)
 	  .attr('width', 22)
 	  .attr('height', 16)
 	  .attr('xlink:href', "/static/icons/seed.gif");
	
	seedRecButton.on("click", function(d) { 
		var url="/nd/make_seed_recommendations/";
		var request = {"pid":d.id}; 
		console.log(d.id); 
		$.post(url, request, function(response) {
  		if (response=='0') {
    			var comment = "Making recommendations based on seed!"; 
    			showRibbon(commend, goodColor, "#mainPanel"); 
    		} else {
      		var failureMessage;
        	if (response=="2") {
          		failureMessage = "Still working on the last batch of recommendations!";
    		} else {
      		failureMessage = "Could not make recommendations from seed!";
    		}			
    		showRibbon(failureMessage, badColor, "#mainPanel"); 
		};
	   })
    ;})
	d3.event.stopPropagation();  
    })
  
 
      nodeElements.attr("transform", function(d) {
	return "translate(" + d.x + "," + d.y + ")";
    });  
  $("#loadingMessage").remove()  
  function grabLinkIndices(links) {
	var linkIndices = [] 
	for (var i=0; i < links.length; i++) { 
		var newSet = {"source": links[i].source.index, "target": links[i].target.index} 
		linkIndices.push(newSet)
        }
	return linkIndices
  }
  var data = {"nodes": JSON.stringify(nodes), "links": JSON.stringify(grabLinkIndices(links))};  
  $.post('/setup_graph/', data);
}, 10);
 });

