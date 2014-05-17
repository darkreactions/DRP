#!/usr/bin/env bash

###########################  User Notes  ###############################
#sh generate_models.sh <database.csv> <compound guide> <trial_name> <parameter_JSON> <OPTIONS...>
# $1 = <database.csv>, $2 = <compound guide>, $3 = <trial_name>, $4 = <parameter_JSON>, $5+ = <OPTIONS>

#Abort the script on unaccounted errors.
set -e

#############################  Check the Input  ########################
#Show "help" and abort if insufficient arguments are given.
NUM_REQ_ARGS="4" #The number of arguments that are required per use.
if [[ $# -lt $NUM_REQ_ARGS ]]; then
	echo "-- Invalid number of arguments! Please use:"
	echo "	sh generate_models.sh <database.csv> <guide.csv> <trial_name> <parameters.json> <OPTIONS...>"
	echo "-- For more information, please see \"README.txt\"."
	echo "" #New line
	exit
fi

#With no options specified, the entire command will run.
CREATE_MODEL=true
CREATE_NEW_HYP=true
RUN_WEKA=true
RUN_WEKA_HYP=true
ANALYZE_HYP=true
GRAPH_MODELS=true
EXISTING_CG=false

if [[ $# -gt $NUM_REQ_ARGS ]]; then
	for arg in ${@:($NUM_REQ_ARGS+1)}
	do
		orig_arg="$arg"
		arg=$(echo "$arg" | tr "[:upper:]" "[:lower:]") #Make lowercase
		arg=$(echo "$arg" | tr -d "-") #Remove hyphens
		arg=$(echo "$arg" | tr -d "_") #Remove underscores
		
		if [ "$arg" = "nonewmodel" ]; then
			CREATE_MODEL=false
		elif [ "$arg" = "nowekamodel" ]; then
			RUN_WEKA=false
		elif [ "$arg" = "nonewhyp" ]; then
			CREATE_NEW_HYP=false
		elif [ "$arg" = "nowekahyp" ]; then
			RUN_WEKA_HYP=false
		elif [ "$arg" = "noanalyzehyp" ]; then
			ANALYZE_HYP=false
		elif [ "$arg" = "nograph" ]; then 
			GRAPH_MODELS=false
		elif [ "$arg" = "nohtml" ]; then 
			GRAPH_MODELS=false
                elif [ "$arg" = "existingcg" ]; then
                        EXISTING_CG=true
		else
			echo "-- Unknown option: $orig_arg "
			echo "-- See README.txt for supported options."
			exit
		fi
	done
fi



############################  Prepare Directories/ETC.  ################
echo "-- Recording run information..."

#Used to display last run-time and for logging files.
DATE=$(date +%F)
DATE=$(echo $DATE| tr "-" "_")

CURRENT_DIR=$(pwd) 
TRIAL_DIR=$CURRENT_DIR/$3
HTML_DIR=$TRIAL_DIR/html
DOT_DIR=$TRIAL_DIR/dots
LOG_DIR=$TRIAL_DIR/log/$DATE

#Make directories that don't exist.
mkdir -p $DOT_DIR
mkdir -p $TRIAL_DIR/sdf

mkdir -p $HTML_DIR/svgs
mkdir -p $HTML_DIR/urls

mkdir -p $LOG_DIR

#Write the command being run and current time to the core log file. 
echo "Run on: $(date)\n" >> $TRIAL_DIR/log/command_log.txt #Date and time.
CLIENT_IP=$SSH_CLIENT
if [[ ! $CLIENT_IP ]] ; then CLIENT_IP="home" ; fi 
echo "By: $(whoami) \@ $CLIENT_IP\n" >> $TRIAL_DIR/log/command_log.txt #SSH IP if applicable.
echo "Arguments: $@\n\n" >> $TRIAL_DIR/log/command_log.txt #Arguments used.

###############################  Building the Model  ###################
if $CREATE_MODEL ; then
	#Construct the expanded csv file.
	echo "-- Expanding and ARFFing..."
        if $EXISTING_CG ; then         
		python scripts/constructDescriptorTable.py $1 $2 $3 existingCG 2> $LOG_DIR/cDT.txt
	else

		python scripts/constructDescriptorTable.py $1 $2 $3 2> $LOG_DIR/cDT.txt
	fi
	wait
	
	#Construct the ARFF files from the expanded.csv.
	python scripts/expanded2arff.py $3 2> $LOG_DIR/exp2arff.txt
	wait
fi

#Run weka to construct initial models.
if $RUN_WEKA ; then
	echo "-- Feeding the weka..."
	python scripts/exportClasspath.py > config/temp_config_file 2> $LOG_DIR/expClassPath.txt
	export CLASSPATH="$CLASSPATH:$(cat config/temp_config_file)"
	rm config/temp_config_file
	
	
	echo "---- Imagining..."
	java weka.filters.unsupervised.instance.RemovePercentage -P 30 -i $TRIAL_DIR/$3_out.arff -o $TRIAL_DIR/$3_out_train.arff
	java weka.filters.unsupervised.instance.RemovePercentage -P 30 -i $TRIAL_DIR/$3_out.arff -o $TRIAL_DIR/$3_out_test.arff -V

	java weka.filters.unsupervised.instance.RemovePercentage -P 30 -i $TRIAL_DIR/$3_pur.arff -o $TRIAL_DIR/$3_pur_train.arff
	java weka.filters.unsupervised.instance.RemovePercentage -P 30 -i $TRIAL_DIR/$3_pur.arff -o $TRIAL_DIR/$3_pur_test.arff -V
	
	echo "---- Realizing..."
	java weka.classifiers.meta.Bagging -d $TRIAL_DIR/$3_pur_bag.model -W weka.classifiers.trees.J48 -t $TRIAL_DIR/$3_pur_train.arff -c last --  > $TRIAL_DIR/$3_pur.tree
	java weka.classifiers.meta.Bagging -d $TRIAL_DIR/$3_out_bag.model -W weka.classifiers.trees.J48 -t $TRIAL_DIR/$3_out_train.arff -c last --  > $TRIAL_DIR/$3_out.tree
	
	echo "---- Fascinating..."
	java weka.classifiers.trees.J48 -T $TRIAL_DIR/$3_pur_test.arff -l $TRIAL_DIR/$3_pur_bag.model -p 0 -c last > $TRIAL_DIR/$3_pur.out
	java weka.classifiers.trees.J48 -T $TRIAL_DIR/$3_out_test.arff -l $TRIAL_DIR/$3_out_bag.model -p 0 -c last > $TRIAL_DIR/$3_out.out
	wait
	
	#Test the models.
	echo "-- Testing models with science..."
	python scripts/model_performance.py $3 
fi

###############################  Graphing the Space  ###################
if $CREATE_NEW_HYP ; then
	#Remember which parameters were used in the calculations.
	cp $4 $TRIAL_DIR/parameters_used.json
	
	#Generate the hypothetical models based off of parameters_used.json.
	echo "-- Hypothesizing and scientificating..."
	python scripts/generate.py $3 2> $LOG_DIR/generate.txt
	wait
	
	echo "-- Calculating and descriptorating..."
	python scripts/constructDescriptorTable.py $TRIAL_DIR/$3_G_expanded.csv $2 $3 existingCG #2> $LOG_DIR/cDT_existingCG.txt 
	wait
	
	echo "-- Cleaning and ARFFing..."
	python scripts/clean2arff.py $3 2> $LOG_DIR/clean2arff.txt
	wait
fi

if $RUN_WEKA_HYP ; then
	python scripts/exportClasspath.py > config/temp_config_file 2> $LOG_DIR/expClassPath.txt
	export CLASSPATH="$CLASSPATH:$(cat config/temp_config_file)"
	rm config/temp_config_file

	echo "-- Asking the mighty weka..."
	java -Xmx12288m weka.classifiers.trees.J48 -T $TRIAL_DIR/$3_H_pur.arff -l $TRIAL_DIR/$3_pur_bag.model -p 0 -c last > $TRIAL_DIR/$3_H_pur.out
	java -Xmx12288m weka.classifiers.trees.J48 -T $TRIAL_DIR/$3_H_out.arff -l $TRIAL_DIR/$3_out_bag.model -p 0 -c last > $TRIAL_DIR/$3_H_out.out
fi

if $ANALYZE_HYP ; then
	#Analyze the hypothetical model with science.
	echo "-- Analyzing the hypothesized model..."
	python scripts/analyze.py $3 2> $LOG_DIR/analyze_errors.txt
	wait
	
	#Build the SMILES file.
	echo "-- Building SMILES..."
	python scripts/constructDescriptorTable.py $1 $2 $3 --validateOnly  2> $LOG_DIR/cDT_validateOnly.txt
	wait
fi

if $GRAPH_MODELS ; then
#Generate similarities and graph the similarities space.
	echo "-- Finding similarities..."
	python scripts/build_similarities.py $3 2> $LOG_DIR/build_similarities.txt
	wait
	
	echo "-- Tossing them around..."
	python scripts/buildclusters.py $3 2> $LOG_DIR/buildclusters.txt
	wait
	
	#Apply and recolor nodes.
	echo "-- Painting le space..."
	python scripts/covered_nodes.py $1 $3 2> $LOG_DIR/covered_nodes.txt 
	wait
	
	python scripts/color.py $3 2> $LOG_DIR/color.txt
	wait
	
	#Construct the graph (svg) with GraphViz.
	echo "-- Playing with GraphViz..."
	neato -Tsvg $DOT_DIR/edges.dot > $HTML_DIR/big_graph.svg 
	wait
fi

if $CREATE_HTML ; then
	##############################  Constructing the HTML  ############
	#Create "a metric ton of triangles" and lump them in a directory "temp."
	echo "-- Making triangles and URLs..."
	python scripts/graph.py $3 2> $LOG_DIR/graph.txt
	wait 
	
	#Do the same with the URLS.
	python scripts/make_urls.py $3 2> $LOG_DIR/make_urls.txt
	wait
fi

#Report successful run-through.
sed -i '$d' $TRIAL_DIR/log/command_log.txt #Remove the line separating entries. 
echo "Success!\n\n" >> $TRIAL_DIR/log/command_log.txt
echo "-- All done! Please see $LOG_DIR/ for log info."

exit
