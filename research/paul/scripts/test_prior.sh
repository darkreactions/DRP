export CLASSPATH=$CLASSPATH:/home/drp/weka/weka.jar

java weka.classifiers.trees.J48 -T /home/drp/research/tmp/$1_out.arff -l /home/drp/research/prior.model -p 0 -c last > /home/drp/research/tmp/$1_prior_out.out

python /home/drp/research/chemml-research-streamlined/scripts/evaluate_prior.py /home/drp/research/tmp/$1_prior_out.out 
