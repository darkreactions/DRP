import json, os, sys

if __name__ == "__main__":

	CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	PREFIX_DIR = "{0}/{1}/".format(CHEMML_DIR, sys.argv[1])	
	HTML_DIR = "{0}/{1}/html/".format(CHEMML_DIR, sys.argv[1])	
	
	hierarchy_json_path = "{0}/hierarchy.json".format(PREFIX_DIR)
	url_list_path = "{0}/urllist.txt".format(PREFIX_DIR)

	urls = None
	hierarchy = None
	with open(url_list_path) as urls:
		urls = json.load(urls)
	with open(hierarchy_json_path) as hier:
		hierarchy = json.load(hier)
		hierarchy = {k: {quad[0]:quad[1:] for quad in hierarchy[k]} for k in hierarchy}

	main_list = []
	for url in urls:
		main_list.append((url,"urls/{}.html".format(url)))
		with open("{}/urls/{}.html".format(HTML_DIR, url), "w") as pair:
			pair.write("<html><body><table><tr><td>Full Name</td><td>First Reactant</td><td>Second Reactant</td><td>Third Reactant</td><td>Link</td></tr>")
			for graph in urls[url]:
				gn,gl = urls[url][graph].split(",")
				print gn
				print gl
				name = gn.split("?")[1]###
				pair.write("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td><a href='../svgs/{}'>{}</a>".format(gn,hierarchy[url][name][0], hierarchy[url][name][1],hierarchy[url][name][2],gl,gl))
			pair.write("</table></body></html>")
	with open("{}/index.html".format(HTML_DIR),"w") as main:
		main.write("<html><body><table>")
		for url in main_list:
			main.write("<tr><td>{}</td><td><a href='{}'>{}</a></td></tr>".format(url[0],url[1],url[1]))
		main.write("</html></body></table>") 

print "---- make_urls.py: Success!"
