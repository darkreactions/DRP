export CLASSPATH=$CLASSPATH:/home/drp/weka/weka.jar

java weka.classifiers.trees.J48 -T $1 -l $2 -p 0 -c last > $3 

