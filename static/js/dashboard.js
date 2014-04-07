d3.json('/get_stats/', function(data) {

//Taken nearly entirely from: http://nvd3.org/examples/lineWithFocus.html 
nv.addGraph(function() {
  //var chart = nv.models.lineWithFocusChart()
  var chart = nv.models.lineWithFocusChart()
                  .x(function(d) { return d[0] })
                  .y(function(d) { return d[1] })
                  .color(d3.scale.category10().range())
                  ;

  //Set up the axes for the actual graph.
  chart.xAxis
      .tickFormat(function(d) {
         return d3.time.format('%b %d')(new Date(d))
      });

  chart.yAxis
      .tickFormat(d3.format(',.1%'));

  //And for the filtering graph...
  chart.x2Axis
      .tickFormat(function(d) {
         return d3.time.format('%b %d')(new Date(d))
      });
  chart.y2Axis
      .tickFormat(d3.format(',.1%'));


  d3.select('#chart svg')
      .datum(data)
      .transition().duration(500)
      .call(chart);

  nv.utils.windowResize(chart.update);

  $("#loadingMessage").remove();

  return chart;
});

});


