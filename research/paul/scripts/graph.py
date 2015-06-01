import json,subprocess,sys, os

CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
PREFIX_DIR = "{0}/{1}/".format(CHEMML_DIR, sys.argv[1])	
DOT_DIR = "{0}/{1}/dots/".format(CHEMML_DIR, sys.argv[1])	
SVG_DIR = "{0}/{1}/html/svgs/".format(CHEMML_DIR, sys.argv[1])		

hierarchy_json_path = "{0}/hierarchy.json".format(PREFIX_DIR)
analyzed_json_path = "{0}/analyzed.json".format(PREFIX_DIR)
url_list_path = "{0}/urllist.txt".format(PREFIX_DIR)

vals = ["400","310","301","220","211","202","130","121","112","103", "040","031","022","013","004"]
def getcolor(pair):
	if pair == 0 or pair[2] == 0:
		return "ffffff"
	color_num = "%2x" % (256 - float(pair[0])/float(pair[2]) * 256)
	color_num_2 = "%2x" % (float(pair[0])/float(pair[2])*256)
	color_str = "ff%s%s" % (color_num, color_num_2)

	return color_str
def clean(s):
	return filter(str.isalnum, str(s))

def build_triangle(ph, orgname):
	pairs = [("400","310"),("400","301"),("310","301"),("310","220"),("310","211"),("301","211"),("301","202"),
		("220","211"),("220","130"),("220","121"),("211","202"),("211","121"),("211","112"),("202","112"),
		("202","103"),("130","121"),("130","040"),("130","031"),("121","112"),("121","031"),("121","022"),("112","103"),
		("112","022"),("112","013"),("103","013"),("103","004"),("040","031"),("031","022"),("022","013"),("013","004")]
	ret_str = ""
	for p in pairs:
		preface = "p"+ph+orgname
		ret_str += preface+p[0] + " -- " + preface+p[1] + "\n"
	return ret_str

table = []
hierarchy = None
with open(hierarchy_json_path) as hier_handle:
	hierarchy = json.load(hier_handle)
info = None
with open(analyzed_json_path) as file_handle:
	info = json.loads(file_handle.read())
urllist = dict()
prefix = sys.argv[1]
for metal_pair in hierarchy: 
	urllist[metal_pair] = dict()
	for reaction in hierarchy[metal_pair]:
		graph = 'graph %s {\n' % reaction[0] 
		for pH in info[reaction[0]]['pH']:
			graph += build_triangle(pH, clean(reaction[0]))
			for pt in info[reaction[0]]['pH'][pH]:
				graph += '%s [style=filled,color="#%s",label="%s"]\n' % ('p' + pH + clean(reaction) + pt, getcolor(info[reaction[0]]['pH'][pH][pt]), pH + "_" + pt)
		graph += '}'
		name = "%s?%s" % (prefix,reaction[0])
		
	#Create the DOT and SVG files.
		with open(DOT_DIR + name + ".dot", 'w') as out_handle:
			out_handle.write(graph)
		subprocess.check_call("neato -Tsvg {0}/{1}.dot > {2}/{1}.svg".format(DOT_DIR, name, SVG_DIR), shell=True)
		
		urllist[metal_pair][name] = "{0},{0}.svg".format(name)

with open(url_list_path,'w') as url:
	json.dump(urllist, url)

print "---- graph.py: Success!"
