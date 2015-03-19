import json, os, subprocess

if __name__ == "__main__":
	
	CHEMML_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
		
	with open(CHEMML_DIR+"/config/config.json", "r") as CONFIG_FILE:
		weka_path = json.loads(CONFIG_FILE.read())["WEKA_JAR_PATH"]

	print weka_path
