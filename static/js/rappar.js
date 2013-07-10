/*!
 * Rappar: SVG to Raphael Parser
 * 
 * Copyright  2011 Room & Board
 * MIT Licensed
 */

(function(){
	
var aa = {
    cursor: 1,
    cx: 0,
    cy: 0,
    fill: 1,
    "fill-opacity": 0,
    font: 1,
    "font-family": 1,
    "font-size": 0,
    "font-style": 1,
    "font-weight": 1,
    height: 0,
    "letter-spacing": 0,
    opacity: 0,
    r: 0,
    rx: 0,
    ry: 0,
    src: 1,
    stroke: 1,
    // "stroke-dasharray": "",
    "stroke-linecap": 1,
    "stroke-linejoin": 1,
    "stroke-miterlimit": 0,
    "stroke-opacity": 0,
    "stroke-width": 0,
    "text-anchor": 1,
    width: 0,
    x: 0,
    y: 0
};

function rappar(svg) {
    var parser = elemental(),
        items = [],
        groups = [],
        idtops = {},
        text,
        textel,
        grad;
    
    function factory(tag, type) {
        return function (data, attr) {
            var el = {type: type || tag, fill: "#000", stroke: "none"};
            items.push(el);
            for (var at in attr) if (attr.hasOwnProperty(at)) {
                eve("rappar." + tag + ".attr." + at, el, attr[at], at);
            }
            eve("rappar." + tag + ".attrend", el, attr);
        }
    }
    
    eve.on("rappar.polygon.attr.points", function (value) {
        this.path = "M" + value + "z";
        eve.stop();
    });
    eve.on("rappar.image.attr.xlink:href", function (value) {
        this.src = value;
        eve.stop();
    });
    eve.on("rappar.path.attr.d", function (value) {
        this.path = value;
        eve.stop();
    });
    eve.on("rappar.*.attr.transform", function (value) {
        this.transform = parseTransform(value);
        eve.stop();
    });
    eve.on("rappar.*.attr.fill", function (value) {
        var id = value.match(/url\(#([^\)]+)\)/),
            el = this;
        if (id) {
            id = id[1];
            if (idtops[id]) {
                this.fill = idtops[id];
            } else {
                eve.on("rappar.found." + id, function (fill) {
                    el.fill = fill;
                });
            }
            eve.stop();
        }
    });
    eve.on("rappar.line.attrend", function (attr) {
        this.path = "M" + [attr.x1, attr.y1, attr.x2, attr.y2];
    });

    eve.on("rappar.*.attrend", function () {
        var i = groups.length;
        while (i--) if (groups[i] != this) {
            for (var key in groups[i]) {
                if (key == "transform" && this.transform) {
                    this.transform = groups[i].transform + this.transform;
                } else {
                    this[key] = groups[i][key];
                }
            }
        }
    });
    eve.on("rappar.*.attr.style", function (value) {
        var f = function () {
            applyStyle(value, this, aa);
            eve.unbind("rappar.*.attrend", f);
        };
        eve.on("rappar.*.attrend", f);
        eve.stop();
    });
    eve.on("rappar.*.attr.*", function (value, name) {
        if (name in aa) {
            this[name] = aa[name] ? value : parseFloat(value);
        }
    });
    eve.on("elemental.tag.circle", factory("circle"));
    eve.on("elemental.tag.ellipse", factory("ellipse"));
    eve.on("elemental.tag.polygon", factory("polygon", "path"));
    eve.on("elemental.tag.path", factory("path"));
    eve.on("elemental.tag.line", factory("line", "path"));
    eve.on("elemental.tag.rect", factory("rect"));
    eve.on("elemental.tag.image", factory("image"));
    eve.on("elemental.tag.text", factory("text"));

    eve.on("rappar.text.attrend", function () {
        text = "";
        textel = this;
        this["text-anchor"] = this["text-anchor"] || "start";
    });
    eve.on("elemental.text", function (data, attr, raw) {
        textel && (text += raw);
    });
    eve.on("elemental./tag.text", function () {
        textel.text = text;
        textel = null;
    });
    eve.on("elemental.tag.g", function (data, attr) {
        var el = {};
        groups.push(el);
        for (var at in attr) if (attr.hasOwnProperty(at)) {
            eve("rappar.g.attr." + at, el, attr[at], at);
        }
        eve("rappar.g.attrend", el);
    });
    eve.on("elemental./tag.g", function (data, attr) {
        groups.pop();
    });
    eve.on("elemental.tag.linearGradient", function (data, attr) {
        grad = {
            id: attr.id,
            angle: +(360 + Raphael.angle(attr.x1, attr.y1, attr.x2, attr.y2, attr.x1 + 100, attr.y1)).toFixed(2),
            stops: []
        };
    });
    eve.on("elemental./tag.linearGradient", function () {
        var s = [grad.angle],
            stop;
        for (var i = 0, ii = grad.stops.length; i < ii; i++) {
            stop = grad.stops[i];
            if (i && i != ii - 1) {
                s.push(stop.offset + ":" + stop.color);
            } else {
                s.push(stop.color);
            }
        }
        idtops[grad.id] = s.join("-");
        eve("rappar.found." + grad.id, null, idtops[grad.id]);
    });
    eve.on("elemental.tag.radialGradient", function (data, attr) {
        grad = {
            id: attr.id,
            stops: []
        };
    });
    eve.on("elemental./tag.radialGradient", function () {
        var s = [],
            stop;
        for (var i = 0, ii = grad.stops.length; i < ii; i++) {
            stop = grad.stops[i];
            if (i && i != ii - 1) {
                s.push(stop.offset + ":" + stop.color);
            } else {
                s.push(stop.color);
            }
        }
        idtops[grad.id] = "r" + s.join("-");
        eve("rappar.found." + grad.id, null, idtops[grad.id]);
    });
    eve.on("elemental.tag.stop", function (data, attr) {
        var stop = {};
        if (attr.style) {
            applyStyle(attr.style, stop);
        }
        stop.offset = stop.offset || attr.offset;
        stop.color = stop["stop-color"] || attr["stop-color"];
        stop.opacity = stop["stop-opacity"] || attr["stop-opacity"];
        if (~(stop.offset + "").indexOf("%")) {
            stop.offset = parseFloat(stop.offset);
        } else {
            stop.offset = parseFloat(stop.offset) * 100;
        }
        stop.offset = +stop.offset.toFixed(2);
        grad.stops.push(stop);
    });
    parser(svg);
    parser.end();
    return items;
}

function parseTransform(t) {
    var m = Raphael.matrix(),
        sep = /\s+,?\s*|,\s*/;
    (t + "").replace(/([a-z]+)\(([^)]+)\)(?:\s+,?\s*|,\s*|$)/gi, function (all, command, values) {
        values = values.split(sep);
        switch (command.toLowerCase()) {
            case "translate":
                m.add(1, 0, 0, 1, values[0], values[1]);
            break;
            case "scale":
                m.add(values[0], 0, 0, values[1], 0, 0);
            break;
            case "rotate":
                m.rotate(values[0], values[1], values[2]);
            break;
            case "skewx":
                m.add(1, 0, Math.tan(Raphael.rad(values[0])), 1, 0, 0);
            break;
            case "skewy":
                m.add(1, Math.tan(Raphael.rad(values[0])), 0, 1, 0, 0);
            break;
            case "matrix":
                m.add(values[0], values[1], values[2], values[3], values[4], values[5]);
            break;
        }
    });
    return m.toTransformString();
}

function applyStyle(css, el, aa) {
    var rules = (css + "").split(";"),
        trim = /^\s+|\s+$/g,
        key;
    for (var i = 0, ii = rules.length; i < ii; i++) {
        var pair = rules[i].split(":");
        key = pair[0].replace(trim, "").replace(/[A-Z]/g, function (letter) {
            return "-" + letter.toLowerCase();
        });
        if (!aa || key in aa) {
            el[key] = pair[1].replace(trim, "");
            if (aa && !aa[key]) {
                el[key] = parseFloat(el[key]);
            }
        }
    }
}

window.rappar = rappar;

})();