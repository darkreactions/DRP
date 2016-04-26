SCRIPT="../build_model.py"
DESCRIPTOR_FILE="descs/new_reactionpH_nonzeroVariance.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
COMMENT="Profiling model building. Legacy rxns. New descriptors including cxcalc. Using reaction pH instead of all pH. 1 random split."
FILE_SUFFIX="legRxnNonZeroCompound_NewCXDscReactionpH_1RandSplit"
SPLITTER="SingleSplitter"

python -u -m cProfile -o "SVMbasic_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_basic -d "SVM $COMMENT" &> "SVMbasic_$FILE_SUFFIX.out" &
