echo "$1"
export CLASSPATH=$CLASSPATH:/home/praccugl/weka-3-6-9/weka.jar
java weka.classifiers.meta.Bagging   -d nodesize_out_bag.model_$1 -t withdata_out_train.arff -W weka.classifiers.trees.J48 -c last -- -M $1 > nodesize.tree_$1 


#java weka.classifiers.trees.J48   -d nodesize_out_bag.model_$1 -t nodesize_out_train.arff  -c last -C .5  > nodesize_pur_out.tree


java weka.classifiers.trees.J48 -T withdata_out_test.arff -l nodesize_out_bag.model_$1 -p 0 -c last > nodesize_out.out_$1

python model_performance.py nodesize $1
