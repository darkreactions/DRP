# First arg is filter result buffer from which feature set should be derived
# Second arg is training set
# Third arg is test set
# Fourth arg is name of file to output model to

bash filter.sh $1 $2 ./J48_train.arff
bash filter.sh $1 $3 ./J48_test.arff
bash J48_make_model.sh ./J48_train.arff ./J48_test.arff $4
rm J48_train.arff
rm J48_test.arff