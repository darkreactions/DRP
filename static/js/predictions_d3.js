var svg = d3.select("svg:nth-child(1)");
var graph = d3.select("#graph1");
var edges = graph.selectAll(".edge").selectAll("path");
var nodes = graph.selectAll(".node").selectAll("ellipse");

var svgWidth = parseFloat(svg.attr("width"));
var svgHeight = svg.getBoundingBox().height;
var halfSvgWidth = svgWidth/2;
var halfSvgHeight = svgHeight/2;

svg
 .attr("preserveAspectRation", "none")
 .attr("width", "100%")
 .attr("height", "100%");
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
 .attr("rx", function(d,i){return (i+1)*15})
 .attr("ry", function(d,i){return (i+1)*15});

nodes
 .on("mouseover", function(d) { //Make the nodes fade in on mouse-in.
  d3.select(this)
   .transition()
   .duration(750)
   .style("stroke","purple")
   .style("fill","cyan")
   .style("stroke-width", "20")
  })
 //Make the nodes fade out on mouse-out.
 .on("mouseout", function(d) {
  d3.select(this)
   .transition()
   .duration(750)
   .style("stroke","orange")
   .attr("stroke-width", "20")
 })
 .on("click", function(d) {
  //Zoom in on the circle that was clicked.
  var cx = this.getBBox().x;
  var cy = this.getBBox().y+svgHeight;  
  var screen_radius = 600;
  var newView = Array(cx-screen_radius, cy-screen_radius, cx+screen_radius, cy+screen_radius);
  svg
   .transition()
   .duration(1000)
   .attr("viewBox", newView.join(" ")); 
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

