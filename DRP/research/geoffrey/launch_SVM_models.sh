SCRIPT="../build_model.py"
DESCRIPTOR_FILE="descs/new_nonZeroVariance.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_reactions"
COMMENT="Default PUK parameters. New Descriptors with NonZero Variance. 15 mutual info splits."
FILE_SUFFIX="PUKdefault_legRxnNonZeroCompound_NewDscNonZeroVar"
SPLITTER="MutualInfoSplitter"

python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_basic -d "SVM $COMMENT" &> "SVMbasic_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_BCR -d "BCR SVM $COMMENT" &> "SVMBCR_$FILE_SUFFIX.out" &
