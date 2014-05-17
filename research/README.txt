########################################################################
#############################__ChemML__#################################
########################################################################


###########################__SETTING UP__###############################
In order to create models, the following need to be installed:
-- Python (2.7)
-- RDkit (for Python)
-- Java
-- GraphViz ("neato" specifically)
-- Weka (3-7-9)
-- ChemAxon (JChem)

Additionally, the configuration file (default: config/config.json) must
	know the location of the .../weka.jar file and .../JChem/bin.

*Note that only Ubuntu 12+ and Weka-3-7-9 have been tested.


#########################__USING ChemML__###############################
To create .svg models and an HTML navigation guide, execute:
"generate_models.sh <database.csv> <guide.csv> <trial_name> <parameters.json>"

Parameter Details:
-- <database.csv> should contain names/amounts/units of reactions as 
		well as information on temperature (celsius), pH, slow-cooling,
		reaction reference, etc..
		
-- <guide.csv> should contain legal translations of reactant 
		abbreviations used in <database.csv>, CAS numbers, and reactant
		types (_i_norganic, _o_rganic, _w_ater, _s_olvent, _p_H).
		
-- <trial_name> is the name of the trial to be run. The hypothetical
		models generated will simply append "_H" to <trial_name>.
		
-- <parameters.json> must be a list of lists and is used to describe the
		the hypothetical space to be generated. The order must be: 
		[WATER MASS in grams], TEMPERATURE, TIME, SLOW-COOLING, pH.
		
		*Example <parameters.json> contents: 
		[[9,27,45], [90,110,130,150],[24,36,48], ["yes"], [1,3,5,7]]
		
		##Note that <parameters.json> will be copied in the 
			trial_directory for record-keeping purposes.


Additional options may be specified after the required arguments above:

OPTIONS:
--NoNewModel	-->  Skips the creation of a new base model.
--NoWekaModel	-->  Skips running Weka commands on the base model.
--NoNewHyp		-->  Skips creating a new hypothetical model.
--NoWekaHyp	-->  Skips running Weka commands on the base model.
--NoAnalyzeHyp	-->  Skips the analysis of the hypothetical model.
--NoGraph		-->  Skips the graph and HTML constructions.
--NoHTML		-->  Skips the graph and HTML constructions.
		
		
#######################__DIRECTORY STRUCTURE__##########################
The scripts assume .../trial_name is a sibling of .../scripts/ and will
	dump files into each trial directory.
	
Thus, to work properly, the ChemML_Directory must resemble:
-- ChemML_Directory/
----- scripts/
-------- manage_models.sh (and all dependent scripts)
-------- mols.json
-------- mlConvert.json
----- config/
-------- config.json
----- trial_name*/		(Each trial_directory is generated per call.)
-------- html/		(Not affected by changing location nor copy/pasting)
----------- svgs 
----------- urls 
----------- index.html 
----------- big_graph.svg 
-------- dots/
-------- sdf/
-------- log/
----------- command_log.txt (Details on each run of this trial.)
----------- [YEAR_MONTH_DAY]/ (Error reports recorded on per-day basis.)


Written by Paul Raccuglia and edited by Casey Falk. 
(COPYRIGHT 2013 -- HAVERFORD COLLEGE)
