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


function zoomed() {
  container.attr("transform",
  "translate(" + d3.event.translate +")scale("+ d3.event.scale + ")");
}


function grabLinkIndices(links) {
  var linkIndices = []
  for (var i=0; i < links.length; i++) {
    var newSet = {
                  "source": links[i].source.index,
                  "target": links[i].target.index
                 }
    linkIndices.push(newSet)
  }
  return linkIndices
}


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

  console.log(links);
  // Size of the text boxes that appear upon hovering over individual nodes.
  var textPlacement = nodeTooltips.length*40;
  var boxLength = 130;

  var translateValue1 = width/2.3
  var translateValue2 = height/2.4
  var scaleValue = height/4500


  var firstCluster = d3.selectAll(".circleClusters1");
  var secondCluster = d3.selectAll(".circleClusters2");
  var tinyNodes = d3.selectAll(".nodeElements");
/*DELETE...also delete dragger.js file!!
  //Here I am creating a zoom slider that will control the zoom level for the svg, and will also tie into which level of the cluster hierarchy is displayed
  var dragger = new Dragdealer('slider', {
      horizontal: false,
      vertical: true,
      animationCallback: function(x, y) {
          $('#slider .value').text(Math.round(y * 100));
          var zoomScaleValue = y + 1 + y*1.5 + 0.1*y*20

          $('#graph').css("transform", "scale(" + zoomScaleValue + ")");

          if (y < 0.3) {
            firstCluster.style("visibility", "visible")
            secondCluster.style("visibility", "hidden")
            tinyNodes.style("visibility", "hidden")
          } else if (y < 0.6) {
            firstCluster.style("visibility", "hidden")
            secondCluster.style("visibility", "visible")
            tinyNodes.style("visibility", "hidden")
          } else if (y < 0.85) {
            firstCluster.style("visibility", "hidden")
            secondCluster.style("visibility", "hidden")
            tinyNodes.style("visibility", "visible")
          }
      }
  });
*/

  var zoom = d3.behavior.zoom()
                .scaleExtent([1/10, 2])
                .on("zoom", zoomed)
                .translate([translateValue1, translateValue2]).scale(scaleValue);


  var drag = d3.behavior.drag()
               .origin(function(d) { return d; })
               .on("dragstart", dragstarted)
               .on("drag", dragged)
               .on("dragend", dragended);

  console.log("1");
  //Add a default weight of   1   for each connection.
  for (var i=0; i<links.length; i++) {
    links[i]["weight"]=1;
  }

  var svg = d3.select("#graph").append("svg")
              .attr("width", width)
              .attr("height", height);

  var container = svg.append("g");

  console.log("2");

  container.attr("transform","translate(" + translateValue1 +"," + translateValue2 + ")scale(" + scaleValue +"," + scaleValue +")");


  // Use a timeout to allow the rest of the page to load first.
    // Run the layout a fixed number of times.
    // The ideal number of times scales with graph complexity.
    // Of course, don't run too long—you'll hang the page!
  var maxIterations = parseFloat(100);
  var loadingBar = $("#innerLoadingBar");
  var loadingBarMaxLength = loadingBar.parent().width();


  console.log("3");
  var force = d3.layout.force()
                .nodes(nodes)
                .links(links)
                .charge(-500)
                .friction(0.6)
                .gravity(0.4)
                .linkDistance(100)
                .size([width, height]);
  console.log("3.5");
  console.log(force);
  console.log(links);

  if (!preLoad) {
    console.log("3.6A")
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
    console.log("3.6B")
    force.start();
    console.log("3.6B2")
    force.tick();
    console.log("3.6B3")
    force.stop();
  }

  console.log("4");


  container.selectAll("line")
    .data(links)
    .enter()
    .append("line")
    .attr("x1", function(d) { return d.source.x; })
    .attr("y1", function(d) { return d.source.y; })
    .attr("x2", function(d) { return d.target.x; })
    .attr("y2", function(d) { return d.target.y; })
    .style("stroke-width", 0.06)
    .attr("stroke", "gray")

  console.log("5");

  container.on("mouseover", function() {
    d3.selectAll(".tooltipContainer").remove();
  });

  var nodeElements = container.selectAll(".node")
                              .data(nodes.filter(function(d) {
                                return d.outcome > 0;
                              }))

  .enter().append("g")
  .attr("class", "node");

  // Clusters of all the nodes with the same SINGLE inorganic in common

  var filter =  container.append("filter").attr("id", "blur")
                         .append("feGaussianBlur")
                         .attr("stdDeviation", 0);



  nodeElements.append("circle")
              .attr("class", "circleClusters1")
              .attr("fill", function(d) {
                return (d.color!=undefined) ? d.color : "rgba(0,0,0,0)";
              })
              .attr("opacity", 0.4)
              .attr("r", 200)
              .attr("filter", "url(#blur)");

//Here I am trying to create a larger, encompassing text element/circle element in order to label the individual clusters(for the single inorg clusters, or circleClusters1
 /*   container.selectAll("circle")
      .data(largeLabels)
      .enter().append("circle")
        .attr("class", "largeLabels")
        .attr("cx", function(d) {return d.x/2;})
        .attr("cy", function(d) {return d.y/2;})
        .attr("r", function(d) {return d.r;})
        .attr("fill", function(d) {return d.fill;});
*/
  //Here I am appending circles that represent the clusters of all the node with the same TWO inorganics in common
    nodeElements.append("circle").attr("class", "circleClusters2").attr("fill", function(d) {return (d.color2!=undefined) ? d.color2 : "rgba(0,0,0,0)";}).attr("opacity", 0.4).attr("r", 80);

    nodeElements.append("circle")
    .attr("class", "nodeElements")
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
	  .attr("in", "SourceGraphic")
	  .attr("stdDeviation", 10)
	  .attr("result", "blur");


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


	var imageContainer = seedRecButton.append("g")
	var seedButton = imageContainer.append("svg:image")
	  .attr('x', textPlacement*1.7)
 	  .attr('y', boxLength + 5)
 	  .attr('width', 22)
 	  .attr('height', 16)
 	  .attr('xlink:href', "/static/icons/seed.png");

	seedButtonRect.on("click", function(d) {
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
    seedButtonRect.append("rect").attr("class", "extraRect").append("text")
		.text("Generate seed recommendations")
        .attr("class", "seedRecText")
		.attr('x', "25px")
		.attr('y', boxLength + 15)

	d3.event.stopPropagation();
    })


      nodeElements.attr("transform", function(d) {
	return "translate(" + d.x + "," + d.y + ")";
    });

    // Hide nodes on the initial zoom level.
    // Upon clicking a circleCluster zoom to that cluster and show the nodes again.
    d3.selectAll(".nodeElements").style("visibility", "visible")
    d3.selectAll(".circleClusters1").style("visibility", "hidden")
    d3.selectAll(".circleClusters2").style("visibility", "hidden")

    d3.selectAll(".circleClusters").on("click", function(d) {
            x = this.cx
            y = this.cy
            shiftX = centerX + x
            shiftY = centerY + y
            container.attr("transform", "translate(" + shiftX + "," + shiftY + ")scale(0.2)")
            d3.selectAll(".nodeElements").style("visibility", "visible")
    ;});

  $("#loadingMessage").remove();

  var data = {
              "nodes": JSON.stringify(nodes),
              "links": JSON.stringify(grabLinkIndices(links))
             };
  $.post('/setup_graph/', data);
 });

