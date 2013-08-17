from django.conf import settings
from django.utils import simplejson

from xml.etree import ElementTree as et
import pygraphviz as pgv
from construct_descriptor_table import *

###
def generate_svg(lab_group, step = "start", source=None):
	###Perform any additional calculations.
	construct_entire_descriptor_table(lab_group)
	
	#Load the edges.dot file:
	overallGraph = pgv.AGraph(settings.DYNAMIC_DIR + "/dots/edges_short.dot")
	overallGraph.layout()
	
	###Replacing seems to have no effect on return time.
	result = overallGraph.draw(format="svg").replace(
		"stroke=\"red\"", "class=\"badNode generalNode\"").replace(
		"stroke=\"yellow\"", "class=\"medNode generalNode\"").replace(
		"stroke=\"blue\"", "class=\"goodNode generalNode\"").replace(
		"style=\"filled\"", "")
	
	return result
	
####FUNCTIONAL, BUT NOT INTEGRATED WITH GRAPHVIZ.
#def generate_svg(lab_group, step = "start", source=None):
	##SVG Settings:
	#svgWidth = "3000"
	#svgHeight = "3000"
	
	##Get the "next step" of the SVG generation.
	#if step=="start":
		#file_name = "compound_list"
		#next_step = "compounds"
	#else:
		#file_name = source #Assumes the path already is stripped of invalid characters.
		#if step=="compounds":
			#next_step = "reactions"
		#else: #if step=="reactions":
			#next_step = "details"
		
		
	##Construct the top of the SVG
	#svg = et.Element("svg", id="svgFile", width=svgWidth, height=svgHeight, viewBox = "0 0 {} {}".format(svgWidth, svgHeight), xmlns="http://www.w3.org/2000/svg")
	
	#attributes = {"id": "graph", "class": "graph {}".format(next_step), "transform":"scale(1 1) rotate(0) translate(4 7872)"}
	#graph = et.SubElement(svg, "g", attrib=attributes)
			
	#try:	
		#with open(settings.DYNAMIC_DIR+"/json/{}.json".format(file_name)) as json_nodes:
			#node_list = simplejson.load(json_nodes)
	#except:
		#raise Exception("Could not load graphic!")
		
	#if next_step != "details":
		#for (node_title, node_attributes) in node_list:
			##Construct the SVG "g"roup node.
			#g = et.SubElement(graph, "g", attrib={"class": "node"})
			#et.SubElement(g, "title").text = node_title #Create the title tag and set its content. 
			#et.SubElement(g, "ellipse", attrib=node_attributes) #Create the ellipse visual.
	
	#else:
		##Construct the mole triangle
		#triangle_points = ["300,-7800", "500,-7500", "100,-7500"]
		#triangle_attributes = {"class":"triangleSVG", "points":" ".join(triangle_points)}
		#g = et.SubElement(graph, "g", attrib={"class": "node"})
		#et.SubElement(g, "polygon", attrib=triangle_attributes) #Create the triangle visual.
		
		##Create anything else the json file specifies (eg, text).
		#i = 0 ###triangle_points
		#text_attributes = {"class":"textNode", "text-anchor":"middle"}
		#text_rotations = ["0","-45","45"]
		#for (node_title, i) in zip(node_list, xrange(3)):
			##Construct the SVG "g"roup node.
			#g = et.SubElement(graph, "g", attrib={"class": "node"})
			##Format and construct the text.
			#text_attributes["x"], text_attributes["y"] = triangle_points[i].split(",")
			
			##Rotate each text element differently.
			#total_attributes = text_attributes
			#total_attributes.update({"transform": "rotate({},{})".format(text_rotations[i], triangle_points[i])})
			#et.SubElement(g, "text", attrib=total_attributes).text = node_title
		
	#return et.tostring(svg)
