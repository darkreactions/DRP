function makeGraph(graphURL, svgContainer) {

  function toFixed(value, precision) {
      var power = Math.pow(10, precision || 0);
      return String(Math.round(value * power) / power);
  }

  //Taken nearly entirely from: http://nvd3.org/examples/lineWithFocus.html
  d3.json(graphURL, function(data) {
    console.log(data)

    nv.addGraph(function() {
      //var chart = nv.models.lineWithFocusChart()
      var chart = nv.models.lineChart()
          .x(function(d) { return d[0]; })
          .y(function(d) { return d[1]; })
          .yDomain([0,1])
          .tooltips(true)
          .tooltipContent(function(key, x, y, e, graph) {
              var HTML = "";
              HTML += '<div class="modelTitle">'+data["titles"][x]+'</div>';
              data["lines"] = data["lines"].sort(function(a,b) {
                if (a["key"]===b["key"]) return 0;
                if (a["key"]>b["key"]) return 1;
                if (a["key"]<b["key"]) return -1;
              });
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

              HTML += '</p><div class="modelDescription">'+data["descriptions"][x] + '</div>';

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

      d3.select(svgContainer)
          .datum(data["lines"])
          .transition().duration(500)
          .call(chart);

      nv.utils.windowResize(chart.update);

      $("#loadingMessage").remove();
      document.querySelectorAll(".nv-x.nv-axis g.nv-axisMaxMin text")[1].classList.add("lastVersion");
      document.querySelectorAll(".nv-x.nv-axis g.nv-axisMaxMin text")[1].style.fontWeight = "bold";

      return chart;
    });
  });
}

makeGraph("/get_class_stats/2","#chart2Class svg");
makeGraph("/get_class_stats/4","#chart4Class svg");
