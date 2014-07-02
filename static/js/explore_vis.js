// D3 code for visualization will go here! 

<!DOCTYPE html>
<meta charset="utf-8">
<style>
 
text {
  font: 10px sans-serif;
}
 
line {
  stroke: #000;
  stroke-width: 1.5px;
}
 
circle {
  stroke: #fff;
  stroke-width: 1.5px;
}

.header{
	text-align: center;
	vertical-align:center;
	background-color: #496388;
	color:#FFF;
	padding:3px;}
	
.body{
	text-indent:1.5em;
	text-align: left;
	vertical-align:center;
	margin-left:5px;
	margin-right:5px;}
	
table {
	margin: 0 auto;
}
 
</style>
<body>
<script src="//ajax.googleapis.com/ajax/libs/jquery/1.11.0/jquery.min.js"></script>

<link rel="stylesheet" href="//ajax.googleapis.com/ajax/libs/jqueryui/1.10.4/themes/smoothness/jquery-ui.css" />
<script src="//ajax.googleapis.com/ajax/libs/jqueryui/1.10.4/jquery-ui.min.js"></script>
<script src="http://d3js.org/d3.v3.min.js"></script>
<script>
  $(function() {
    $(document).tooltip();
	for (var i = 1; i < 5; i++) {
		console.log("HERE");
          $("td#color"+i).css({
          	
            "background-color": function(){
            	console.log(i);
            	if (i == 4)	return "#1a9641";
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
  
  
</script>
<script>
var width = 1800,
    height = 2000;
 


d3.json("data", function(graph) {
    var n = 100,
    nodes = graph.nodes,
    links = graph.links;
	
	var force = d3.layout.force()
    .nodes(nodes)
    .links(links)
    .charge(-350)
    .friction(0.6)
    .gravity(0.4)
    .linkDistance(4)
    .size([width, height]);
 
var svg = d3.select("body").append("svg")
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
 
  
</script>
<div class="header"><h2>Dark Reactions Network</h2></div>
        <div class="body">
        <p>	Each node is linked to its 5 nearest neighbors in the graph. The PageRank of each node is then calculated, and displayed the graph as a force directed graph. The colors indicate the outcome of the particular experiment. 
        </p>
        <p>This network visually groups experiments that are similar. The PageRank alone is not particularly  useful, as it only indicates whether the experiment is similar to many other experiments. However, by looking at nodes with small PageRank and high outcome, one could seek potentially useful new avenues of exploration. Conversely, looking at dense groups of low outcomes could give a scientist an idea of what kinds of experiments are not worth reproducing.
        </p>
        </div>
	</div>
	
        <table>
            <td>Outcome 1</td>
            <td>Outcome 2</td>
            <td>Outcome 3</td>
            <td>Outcome 4</td>
          </tr>
          <tr>
            <td id="color1"></td>
            <td id="color2"></td>
            <td id="color3"></td>
            <td id="color4"></td>
          </tr>
        </table>
</body>
    	

