d3.json('/get_stats/', function(data) {

//Taken nearly entirely from: http://nvd3.org/examples/lineWithFocus.html
nv.addGraph(function() {
  //var chart = nv.models.lineWithFocusChart()
  var chart = nv.models.lineChart()
      .x(function(d) { return d[0]; })
      .y(function(d) { return d[1]; })
      .yDomain([0,1])
      .tooltips(true)
      .tooltipContent(function(key, x, y, e, graph) {
          var HTML = "";
          for (i=0; i< data["lines"].length; i++) {
              var cur_key = data["lines"][i]["key"];
              var values = data["lines"][i]["values"];
              var versionSpecific = values.filter(function(tup){
                return tup[0]==x;
              });
              if (versionSpecific.length) {
                var val = versionSpecific[0][1];
                if (val<=1.0 && val>=0.0) {
                  HTML += cur_key + ': ' + toFixed(val*100, 1) + '%<br/>';
                } else {
                  HTML += cur_key + ': ' + val + '<br/>';
                }
              }
            };

          HTML += '</p><p> Model notes:<br />' + data["model_labels"][x] + '</p>';

          return HTML
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
  document.querySelectorAll(".nv-x.nv-axis g.nv-axisMaxMin text")[1].classList.add("lastVersion");
  document.querySelectorAll(".nv-x.nv-axis g.nv-axisMaxMin text")[1].style.fontWeight = "bold";

/* TODO: Rewrite this mess in straight D3 */

/* TODO: An array with the appropriate descriptions is in data["model_labels"]
However, NVD3 (the library of nifty charts) doesn't seem to play nicely with
custom tooltips...
*/

  return chart;
});

});

function toFixed(value, precision) {
    var power = Math.pow(10, precision || 0);
    return String(Math.round(value * power) / power);
}
