var svg = d3.select("svg:nth-child(1)");
var graph = d3.select("#graph1");
var edges = graph.selectAll(".edge").selectAll("path");
var nodes = graph.selectAll(".node").selectAll("ellipse");

var svgWidth = parseFloat(svg.attr("width"));
var svgHeight = parseFloat(svg.attr("height"));
var halfSvgWidth = svgWidth/2;
var halfSvgHeight = svgHeight/2;

svg
 .attr("preserveAspectRation", "none")
 .attr("width", "100%")
 .attr("height", "100%")
 /* Mouse Navigation? Not stable/functional.
 .on("click", function(d) {
  var newLocation = d3.mouse(this);
  alert(newLocation);
  var newViewBox = Array(
   newLocation[0]-halfSvgWidth,
   newLocation[1]-halfSvgHeight,
   newLocation[0]+halfSvgWidth,
   newLocation[1]+halfSvgHeight
   );
  alert(newViewBox);
  svg
   .transition()
   .duration(1000)
   .attr("viewBox", newViewBox.join(" "));
});
*/


nodes
 .transition()
 .duration(1000)
 .attr("rx", function(d,i){return (i+1)*5})
 .attr("ry", function(d,i){return (i+1)*5})

  //Make the nodes fade in on mouse-in.
nodes
 .on("mouseover", function(d) {
 d3.select(this)
  .transition()
  .duration(750)
  .style("stroke","green")
  .style("fill","green");
})
 //Make the nodes fade out on mouse-out.
 .on("mouseout", function(d) {
 d3.select(this)
  .transition()
  .duration(750)
  .attr("rx", 0)
  .attr("ry",0)
  .remove();
})
 .on("click", function(d) {
 d3.select(this)
  .transition()
  .duration(400)
  .style("fill", "cyan")
  .style("stroke", "navy")
  .attr("rx", 100)
  .attr("ry", 100);
});

edges
 .on("mouseover", function(d) {
 d3.select(this)
  .transition()
  .duration(500)
  .style("fill", "black")
  .style("stroke", "black")
  .style("stroke-width", "20");
})
 //.on("mouseout", function(d) {
 //d3.select(this)
  //.transition()
  //.duration(500)
  //.style("fill", "grey")
  //.style("stroke-width", "0");
//});
;

