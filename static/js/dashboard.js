d3.json('/get_stats/', function(data) {

//Taken nearly entirely from: http://nvd3.org/examples/lineWithFocus.html 
nv.addGraph(function() {
  //var chart = nv.models.lineWithFocusChart()
  var chart = nv.models.lineChart()
                  .x(function(d) { return d[0] })
                  .y(function(d) { return d[1] })
                  .yDomain([0,1])
                  //.useInteractiveGuideline(true) // Can't use this and customized tooltips simultaneously
                  .tooltips(true)
                  .tooltipContent(function(key, x, y, e, graph) {
                      return '<h3>' +  key + '</h3>' + '<p>' + y + ' at ' + x  + '</p>';
                  })
                  .showXAxis(true)
                  ;

  //Set up the axes for the actual graph.
  chart.xAxis
      .axisLabel("Version")
      .tickFormat(d3.format(",r"));

  chart.yAxis
      .tickFormat(d3.format(',.1%'));

  d3.select('#chart svg')
      .datum(data["lines"])
      .transition().duration(500)
      .call(chart);

  nv.utils.windowResize(chart.update);

  $("#loadingMessage").remove();

/* TODO: An array with the appropriate descriptions is in data["model_labels"]
However, NVD3 (the library of nifty charts) doesn't seem to play nicely with
custom tooltips...
*/

  return chart;
});

});


