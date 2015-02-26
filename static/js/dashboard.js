function makeGraph(graphURL, svgContainer) {

  function toFixed(value, precision) {
      var power = Math.pow(10, precision || 0);
      return String(Math.round(value * power) / power);
  }

  function formatConfusionTable(jsonTable) {
    var HTML = '<div class="modelConfusionTable">';

    HTML += "<table>";
    HTML += "<tr><td></td><td colspan=5><span class='confusionTableLabel'>Guesses</span></td></tr>";
    HTML += "<tr><td rowspan=6><span class='confusionTableLabel'>Actual</span></td></tr>";

    for (var i=0; i < jsonTable.length; i++) {

      HTML += "<tr>";
      for (var j=0; j < jsonTable.length; j++) {
        HTML += "<td>"+jsonTable[i][j]+"</td>";
      }
      HTML += "</tr>";

    }

    HTML += "</table>";

    HTML += "</div>";
    return HTML;
  }

  //Taken nearly entirely from: http://nvd3.org/examples/lineWithFocus.html
  d3.json(graphURL, function(data) {


    nv.addGraph(function() {
      var chart = nv.models.lineChart()
          .x(function(d) { return d[0]; })
          .y(function(d) { return d[1]; })
          .tooltips(true)
          .tooltipContent(function(key, x, y, e, graph) {
              var HTML = "";

              HTML += '<div class="modelTitle">'+data["titles"][x]+'</div>';
              HTML += '<div class="modelDescription">'+data["descriptions"][x] + '</div>';

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

              HTML += formatConfusionTable(data["confusionTables"][x]);

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

      // Remove the "Test Size" line from the lines to graph.
      var lines = data["lines"].filter(function(obj) {
        return obj["key"]!=="Test Size" && obj["key"]!=="Train Size";
      });


      d3.select(svgContainer)
          .datum(lines)
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

var query2 =  $("#chart2Class").attr("query");
var query4 =  $("#chart4Class").attr("query");

makeGraph("/get_class_stats/2-test"+query2,"#chart2Class svg");
makeGraph("/get_class_stats/4-test"+query4,"#chart4Class svg");
