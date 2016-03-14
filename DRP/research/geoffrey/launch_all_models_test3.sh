SCRIPT="build_model.py"
DESCRIPTOR_FILE="descs/small.dsc"
OUTCOME_DESCRIPTOR="boolean_crystallisation_outcome"
RXN_SET_NAME="hello_testing_datasets_0"
COMMENT="Test3"
FILE_SUFFIX="Test3"
SPLITTER="KFoldSplitter"

python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_basic -d "SVM $COMMENT" &> "SVMbasic_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_BCR -d "$BCR SVM COMMENT" &> "SVMBCR_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt J48 -d "J48 $COMMENT" &> "J48_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt KNN -d "KNN $COMMENT" &> "KNN_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt NaiveBayes -d "NB $COMMENT" &> "NB_$FILE_SUFFIX.out" &
