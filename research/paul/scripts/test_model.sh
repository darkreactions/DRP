export CLASSPATH=$CLASSPATH:/home/drp/weka/weka.jar

java weka.classifiers.trees.J48 -T /home/drp/research/tmp/$1_out.arff -l /home/drp/research/2.3.14_svm_true.model -p 0 -c last > /home/drp/research/tmp/$1_out.out
python /home/drp/research/chemml-research-streamlined/scripts/evaluate_outcomes.py /home/drp/research/tmp/$1_out.out 


