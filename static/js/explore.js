// D3 code for visualization will go here!

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
  }
  );


var width = 1800,
    height = 2000;

d3.json("/get_graph/", function(graph) {
    var n = 100;
    var nodes = graph.nodes;
    var links = graph.links;

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
                 .charge(-350)
                 .friction(0.6)
                 .gravity(0.4)
                 .linkDistance(4)
                 .size([width, height]);

var svg = d3.select("#graph").append("svg")
    .attr("width", width)
    .attr("height", height);


var loading = svg.append("text")
    .attr("x", width / 2)
    .attr("y", height / 2)
    .attr("dy", ".35em")
    .style("text-anchor", "middle")
    .text("Simulating. One moment please…");

// Use a timeout to allow the rest of the page to load first.
setTimeout(function() {

  // Run the layout a fixed number of times.
  // The ideal number of times scales with graph complexity.
  // Of course, don't run too long—you'll hang the page!
  force.start();
  for (var i = n*n; i > 0; --i) force.tick();
  force.stop();

  svg.selectAll("line")
    .data(links)
  .enter().append("line")
    .attr("x1", function(d) { return d.source.x; })
    .attr("y1", function(d) { return d.source.y; })
    .attr("x2", function(d) { return d.target.x; })
    .attr("y2", function(d) { return d.target.y; })
    .style("stroke-width", 0.06)
     .attr("stroke", "gray");

  svg.selectAll("circle")
    .data(nodes.filter(function(d) { return d.outcome > 0;}))
    .enter().append("circle")
    .attr("cx", function(d) { return d.x; })
    .attr("cy", function(d) { return d.y; })
    .attr("r", function(d) { var size = Math.abs(Math.log(d.pagerank))/3 + d.pagerank*4500;
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
    .call(force.drag)
    .append("title").text(function(d) {return "Title: " + d.name + "\nPagerank: " + (d.pagerank).toFixed(5) + "\nPurity: " + d.purity + "\nOutcome: " + d.outcome + "\ninorg1: " + d.inorg1 + "\ninorg2: " + d.inorg2 + "\norg1: " + d.org1;});


  loading.remove();
}, 10);
 });


