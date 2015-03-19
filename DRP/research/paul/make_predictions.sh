export CLASSPATH=$CLASSPATH:/home/drp/weka/weka.jar

java weka.classifiers.trees.J48 -T $1 -l $2 -c last -p 0 > $3
