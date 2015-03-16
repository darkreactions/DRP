export CLASSPATH=$CLASSPATH:/home/drp/weka/weka.jar

java weka.classifiers.trees.J48 -T /home/drp/web/darkreactions.haverford.edu/app/DRP/tmp/$1_out.arff -l /home/drp/research/2.3.14_svm_true.model -p 0 -c last > /home/drp/web/darkreactions.haverford.edu/app/DRP/tmp/$1.out


