// D3 code for visualization will go here!

var width = parseInt(d3.select("#graph").style("width"), 10);
var height = parseInt(d3.select("#graph").style("height"), 10); 
console.log(width + "," + height) 

d3.json("/get_graph/", function(graph) {
    
    var nodeTooltips = [["Reference Number", "ref"],
			["Pagerank", "pagerank"],
			["Purity", "purity"],
			["Outcome", "outcome"],
			["Inorg1", "inorg1"],
			["Inorg2", "inorg2"],
			["Org1", "org1"] 
			]; 
    var n = 100;
    var nodes = graph.nodes;
    var links = graph.links;

    var zoom = d3.behavior.zoom()
    .scaleExtent([1/10, 2])
    .on("zoom", zoomed);

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

var container = svg.append("g")

container.append("rect")
	.attr("width", width)
	.attr("height", height)
	.attr("fill", "white") 
	.on("click", function() { console.log("Clicked the container"); d3.selectAll(".tooltipContainer").remove();}); 


 
function zoomed() {
  container.attr("transform", 
  "translate(" + d3.event.translate +")scale("+ d3.event.scale + ")");
}

// Use a timeout to allow the rest of the page to load first.
setTimeout(function() {

  // Run the layout a fixed number of times.
  // The ideal number of times scales with graph complexity.
  // Of course, don't run too long—you'll hang the page!
 force.on("tick", function() {
  nodes[0].x = width / 2;
  nodes[0].y = height / 2; 
}); 

  force.start();
  for (var i = n*n; i > 0; --i) force.tick();
  force.stop();

  container.selectAll("line")
    .data(links)
  .enter().append("line")
    .attr("x1", function(d) { return d.source.x; })
    .attr("y1", function(d) { return d.source.y; })
    .attr("x2", function(d) { return d.target.x; })
    .attr("y2", function(d) { return d.target.y; })
    .style("stroke-width", 0.06)
     .attr("stroke", "gray");

  var nodeElements = container.selectAll(".node")
    .data(nodes.filter(function(d) { return d.outcome > 0;}))

    .enter().append("g")
    .attr("class", "node"); 
 
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
    .on("mousedown", function() {
	d3.event.stopPropagation()})

    .on("click", function() {
	d3.selectAll(".tooltipContainer").remove(); 
	var thisGroup = this.parentNode;
	this.parentNode.parentNode.appendChild(thisGroup); 
	
	var currentCircle = d3.select(thisGroup);  
	console.log("Made it to mouseover")  
	
	var textbox = currentCircle.append("g")
	  .attr("class", "tooltipContainer") 
        
	textbox.append("rect")
	  .attr("class", "tooltipBackground") 
	  .attr("width", 300)
	  .attr("height", nodeTooltips.length*35);
	
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
		   .attr("dy", "2em") 
		}
	textbox.append("text").text("Click to generate seed recs: ") 
		.attr("id", "seedRecText") 
		.attr('x', "20px") 
		.attr('y', nodeTooltips.length*32) 
 
	var imageContainer = textbox.append("g")
	var seedButton = imageContainer.append("svg:image")
 	  .attr("id", "seedImage") 
	  .attr('x', nodeTooltips.length*30) 
 	  .attr('y', nodeTooltips.length*30)
 	  .attr('width', 20)
 	  .attr('height', 20)
 	  .attr('xlink:href', "/static/icons/seed.gif")
	  .on("click", function(d) { 
		var url="/make_seed_recommendations/";
		var request = {"pid":d.id}
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
	   });}) 
	d3.event.stopPropagation();  
    }) 
    nodeElements.attr("transform", function(d) {
	return "translate(" + d.x + "," + d.y + ")";
    });  
  $("#loadingMessage").remove()   
}, 10);
 });

