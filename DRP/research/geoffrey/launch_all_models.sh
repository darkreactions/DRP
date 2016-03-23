SCRIPT="../build_model.py"
DESCRIPTOR_FILE="descs/new_full.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
COMMENT="Legacy rxns. New descriptors NOT including cxcalc. 1 mutual info split."
FILE_SUFFIX="legRxnNonZeroCompound_NewNOCXDsc_1MISplit"
SPLITTER="MutualInfoSplitter"

python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_basic -d "SVM $COMMENT" &> "SVMbasic_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_BCR -d "BCR SVM $COMMENT" &> "SVMBCR_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt J48 -d "J48 $COMMENT" &> "J48_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt KNN -d "KNN $COMMENT" &> "KNN_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt NaiveBayes -d "NB $COMMENT" &> "NB_$FILE_SUFFIX.out" &