export CLASSPATH=$CLASSPATH:/home/drp/weka/weka.jar


java weka.classifiers.functions.SMO -d $1 -t $2 -c last -K "weka.classifiers.functions.supportVector.Puk -O 0.5 -S 7" -p 0
