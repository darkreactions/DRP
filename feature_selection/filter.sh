featureset=`bash parse_feature_set.sh $1`
java -cp /home/drp/weka/weka.jar weka.filters.unsupervised.attribute.Remove -V -R $featureset -i $2 -o $3