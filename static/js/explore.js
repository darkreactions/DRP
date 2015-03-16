$(document).ready(function() {

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


console.log(parseInt(d3.select("#graph").style("width")))
var width = parseInt(d3.select("#graph").style("width"), 10);
var height = parseInt(d3.select("#graph").style("height"), 10);
d3.json("/get_graph/", function(graph) {
var centerX = width/2
var centerY = height/2
  var nodeTooltips = [["Purity", "purity"],
      		["Outcome", "outcome"],
      		["Reference Number", "ref"],
      		["Inorg1", "inorg1"],
      		];

  var nodes = graph.nodes;
  //var clusterNodes = graph.clusterNodes;
  var preLoad = graph.skipTicks === "True";
  var links = graph.links;
  var clusters1 = graph.clusters1;
  var clusters2 = graph.clusters2; 
  // Size of the text boxes that appear upon hovering over individual nodes.
  var textPlacement = nodeTooltips.length*40;
  var boxLength = 130;

  var zoomLevel = "circles1"
  var translateValue1 = width/3.2
  var translateValue2 = height/2.3 //this controls up and down
  var scaleValue = height/2700 
  console.log(translateValue2)
  var firstCluster = d3.selectAll(".circleClusters1");
  var secondCluster = d3.selectAll(".circleClusters2");
  var tinyNodes = d3.selectAll(".nodeElements");

  console.log("1");
  //Add a default weight of   1   for each connection.
  for (var i=0; i<links.length; i++) {
    links[i]["weight"]=1;
  }

  var svg = d3.select("#graph").append("svg")
              .attr("width", width)
              .attr("height", height)

  var container = svg.append("g") 

  console.log("2");

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

 
  var zoomInitial = d3.behavior.zoom()
                .scaleExtent([1/20, 2])
                .on("zoom", zoomed)
                .translate([translateValue1, translateValue2]).scale(scaleValue);

  var zoom = d3.behavior.zoom()
        .translate([0,0])
        .scale(1)
        .scaleExtent([1/2, 8])
        .on("zoom", zoomed);

  var drag = d3.behavior.drag()
               .origin(function(d) { return d; })
               .on("dragstart", dragstarted)
               .on("drag", dragged)
               .on("dragend", dragended);


  svg.call(zoom).call(zoom.event); 

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

  container.on("mouseover", function() {
    d3.selectAll(".tooltipContainer").remove();
  });

  var tip = d3.tip()
      .attr("class", "d3-tip")
      .offset([-10,0])
      .html(function(d) {
        return d.label;
      });

svg.call(tip);

    // Clusters of all the nodes with the same SINGLE inorganic in common
  function createCircleClusters1() {
        var circleClusters1 = container.selectAll("circle")
              .data(clusters1)
              .enter().append("circle") 
              .attr("class", "circleClusters1")
              .attr("fill", function(d) {
                return (d.color!=undefined) ? d.color : "rgba(0,0,0,0)";
              })
              .attr("opacity", 0.4)
              .attr("cx", function(d) { return d.x;})
              .attr("cy", function(d) { return d.y;}) 
              .attr("r", 40) 
              .on("mouseover", tip.show)
              .on("mouseout", tip.hide) 
              .on("click", clicked); 
  }; 


  //Clusters of all nodes with both inorganics in common  
  function createCircleClusters2() {
              var circleClusters2 = container.selectAll("circle")
              .data(clusters2)
              .enter().append("circle") 
              .attr("class", "circleClusters2")
              .attr("fill", function(d) {
                return (d.color!=undefined) ? d.color : "rgba(0,0,0,0)";
              })
              circleClusters2
              .attr("opacity", 0.4)
              .attr("cx", function(d) { return d.x;})
              .attr("cy", function(d) { return d.y;}) 
              .attr("r", 200) 
              .on("mouseover", tip.show)
              .on("mouseout", tip.hide) 
              .on("click", clicked);
   }; 

//Nodes representing each individual reactions 
function createNodeElements() {
  container.selectAll("line")
     .data(links)
     .enter()
     .append("line")
     .attr("class", "lines") 
     .attr("x1", function(d) { return d.source.x; })
     .attr("y1", function(d) { return d.source.y; })
     .attr("x2", function(d) { return d.target.x; })
     .attr("y2", function(d) { return d.target.y; })
     .style("stroke-width", 0.06)
     .attr("stroke", "gray") 
  
   var baseNodes = container.selectAll("g")
    .data(nodes.filter(function(d) {
            return d.outcome > 0;
    }))
    .enter().append("g")
    .attr("class", "node");

    baseNodes.append("circle")
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

	var seedButtonRect = seedRecButton.append("rect")
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

    baseNodes.attr("transform", function(d) {
	return "translate(" + d.x + "," + d.y + ")";
    });
}

function reset() {
    console.log("reset")
    zoomLevel = "first"
    container.attr("transform","translate(" + translateValue1 +"," + translateValue2 + ")scale(" + scaleValue +"," + scaleValue +")");
    $(".links").remove() 
    $(".nodeElements").remove() 
    $(".circleClusters1").remove()
    $(".label1").remove()
    $(".label2").remove()
    $(".circleClusters2").remove()  
    createCircleClusters1();
    }; 


  
function clicked(d) {
        if (zoomLevel == "first") {
          zoomLevel = "second" 
          var dx = Math.round(d3.select(this).attr("r")) * 2,
              dy = dx,
              x = Math.round(d3.select(this).attr("cx")),
              y = Math.round(d3.select(this).attr("cy")),
              scale = .9 / Math.max( dx / width, dy / height),
              translate = [width /2 - scale* x, height/2 - scale*y];
          $(".circleClusters1").remove()
          $(".label1").remove()
          svg.transition()
            .duration(750)
            .call(zoom.translate(translate).scale(scale).event);
          createCircleClusters2();
        } else if (zoomLevel === "second") {
          zoomLevel = "third" 
          console.log(zoomLevel) 
          var dx = Math.round(d3.select(this).attr("r")) * 2,
              dy = dx,
              x = Math.round(d3.select(this).attr("cx")),
              y = Math.round(d3.select(this).attr("cy")),
              scale = .9 / Math.max( dx / width, dy / height),
              translate = [width /2 - scale* x, height/2 - scale*y];
          $(".circleClusters2").remove()
          $(".label2").remove()
          svg.transition()
            .duration(750)
            .call(zoom.translate(translate).scale(scale).event);
         createNodeElements(); 
         } else {
         return reset() 
         }
}


       
  $("#reset").on("click", reset);

  $("#loadingMessage").remove();
  container.attr("transform","translate(" + translateValue1 +"," + translateValue2 + ")scale(" + scaleValue +"," + scaleValue +")"); 
  //createCircleClusters2();  
  //createLabel2s();
  //createNodeElements();
  createCircleClusters1(); 
  var data = {
              "nodes": JSON.stringify(nodes),
              "links": JSON.stringify(grabLinkIndices(links))
             };
  $.post('/setup_graph/', data);
 });
}); 
