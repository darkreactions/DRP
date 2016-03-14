SCRIPT="../build_model.py"
DESCRIPTOR_FILE="../descs/small.dsc"
OUTCOME_DESCRIPTOR="boolean_crystallisation_outcome"
RXN_SET_NAME="hello_testing_datasets_0"
COMMENT="TEST"
FILE_SUFFIX="Test"

python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s MutualInfoSplitter -v -mt SVM_PUK_basic -d "$COMMENT SVM, 15 mutualinfo splits" &> "SVMbasic_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s MutualInfoSplitter -v -mt SVM_PUK_BCR -d "$COMMENT BCR SVM, 15 mutualinfo splits" &> "SVMBCR_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s MutualInfoSplitter -v -mt J48 -d "$COMMENT J48, 15 mutualinfo splits" &> "J48_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s MutualInfoSplitter -v -mt KNN -d "$COMMENT KNN, 15 mutualinfo splits" &> "KNN_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s MutualInfoSplitter -v -mt NaiveBayes -d "$COMMENT NB, 15 mutualinfo splits" &> "NB_$FILE_SUFFIX.out" &
