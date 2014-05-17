echo "$1"
export CLASSPATH=$CLASSPATH:/home/praccugl/weka-3-6-9/weka.jar
java weka.classifiers.meta.Bagging   -d conf_out_bag.model_$1 -t withdata_out_train.arff -W weka.classifiers.trees.J48 -c last -- -C $1 > nodesize.tree_$1 
java weka.classifiers.trees.J48 -T withdata_out_test.arff -l conf_out_bag.model_$1 -p 0 -c last > conf_out.out_$1

python model_performance.py conf $1
