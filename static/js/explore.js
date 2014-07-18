  $(function() {
    $(document).tooltip();
      for (var i = 1; i < 5; i++) {
        console.log("HERE");
        $("td#color"+i).css({
            "background-color": function(){
              console.log(i);
              if (i == 4)  return "#1a9641";
              else if (i == 3) return "#a6d96a";
              else if (i == 2) return "#fdae61";
              else if (i == 1) return "#d7191c";
              },
            "width": "20px",
            "height": "20px"
        })
     }
  });


  var width = parseInt(d3.select("#graph").style("width"), 10);
  var height = parseInt(d3.select("#graph").style("height"), 10); 

d3.json("/get_graph/", function(graph) {
    var n = 100;
    var nodes = graph.nodes;
    var links = graph.links;
    var dids = graph.dids;  
   
    var zoom = d3.behavior.zoom()
    .scaleExtent([1/10, 2])
    .on("zoom", zoomed);

    var drag = d3.behavior.drag()
      .origin(function(d) { return d; })
      .on("dragstart", dragstarted)
      .on("drag", dragged)
      .on("dragend", dragended);

/*Function to bring up 'seed-rec' icon on mouseover, make it disappear on mouseout
function mouseover(d) { 
  d3.select(this)
  .classed("dataSpecificButtonContainer", true); 
}

*/
//styling for toolTips
  function prettifiedTooltip(tooltips, i, d){
    var title = tooltips[i][1] + ": ";
    var rawText = d[tooltips[i][0]];
    var prettyText = rawText.split("+").join(", ");
    return title + prettyText;
  }
 
//Needed for zooming and dragging (http://bl.ocks.org/mbostock/6123708).
  function dragstarted(d) {
    d3.event.sourceEvent.stopPropagation();
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

  var container = svg.append("g"); 

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

  container.selectAll("circle")
    .data(nodes.filter(function(d) { return d.outcome > 0;}))
    .enter().append("circle")
    .attr("cx", function(d) { return d.x; })
    .attr("cy", function(d) { return d.y; })
    .attr("r", function(d) { var size = Math.abs(Math.log(d.pagerank))/3 + d.pagerank*1000;
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
    .on("mouseover", function() {
        var g = this.parentNode;
        this.parentNode.parentNode.appendChild(g)
        d3.select(g).select("circle");
        
	var textbox = d3.select(g).append("g")
            .attr("class", "tooltipContainer")

	textbox.append("rect")
        .attr("x", -100)
 	.attr("y", -20*tooltips.length-15)
	.attr("width", 400)
	.attr("height", 20*tooltips.length+10)
	
	for (var i=0; i < tooltips.length; i++){
		textbox.append("text") 
		.attr("dx", -80)
		.attr("dy", (-20*tooltips.length+20)*i+5)
		.attr("text-anchor", "left")
		.attr("class", "tooltip")
		.text(function(d) {return "Reference Number: " + d.ref + "\nPagerank: " + (d.pagerank).toFixed(5) + "\nPurity: " + d.purity + "\nOutcome: " + d.outcome + "\ninorg1: " + d.inorg1 + "\ninorg2: " + d.inorg2 + "\norg1: " + d.org1;
		}); 
	} 
	})
	.on("mouseout", function() {
		d3.select(".tooltipContianer").remove(); 

  $("#loadingMessage").remove()   
}, 10);
});

/*Put this on container that pops up on click event: 
function mouseover() {

 var url="/make_seed_recommendations/";
 var request = {"did":did} 

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
*/ 
