from xml.etree import ElementTree as et
from django.conf import settings
from django.utils import simplejson

def generate_svg(lab_group):
	svgWidth = "3000"
	svgHeight = "3000"
	#Construct the top of the SVG
	svg = et.Element("svg", id="svgFile", width=svgWidth, height=svgHeight, viewBox = "0 0 {} {}".format(svgWidth, svgHeight), xmlns="http://www.w3.org/2000/svg")
	
	attributes = {"id": "graph1", "class": "graph", "transform":"scale(1 1) rotate(0) translate(4 7872)"}
	graph = et.SubElement(svg, "g", attrib=attributes)
			
	with open(settings.DYNAMIC_DIR+"/json/node_list.json") as json_nodes:
		node_list = simplejson.load(json_nodes)
	
	for (node_title, node_attributes) in node_list:
		#Construct the SVG "g"roup node.
		g = et.SubElement(graph, "g", attrib={"class": "node"})
		et.SubElement(g, "title").text = node_title #Create the title tag and set its content. 
		et.SubElement(g, "ellipse", attrib=node_attributes) #Create the ellipse visual.
		
	return et.tostring(svg)
