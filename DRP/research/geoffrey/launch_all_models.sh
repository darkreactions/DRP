#!/bin/bash
SCRIPT="../build_model.py"
DESCRIPTOR_FILE="descs/legacy_nonZeroVariance.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
COMMENT="Legacy rxns. Legacy descriptors NonzeroVariance missing != -1. 15 mutual info split."
FILE_SUFFIX="legRxnNonZeroCompound_LegacyDscNonzeroVarMissingNone_15MISplit"
SPLITTER="MutualInfoSplitter"

python -u -m cProfile -o "SVMbasic_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_basic -d "SVM Default PUK $COMMENT" &> "SVMbasic_$FILE_SUFFIX.out" &
python -u -m cProfile -o "SVMBCR_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_BCR -d "BCR SVM Default PUK $COMMENT" &> "SVMBCR_$FILE_SUFFIX.out" &
python -u -m cProfile -o "J48_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt J48 -d "J48 $COMMENT" &> "J48_$FILE_SUFFIX.out" &
python -u -m cProfile -o "KNN_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt KNN -d "KNN $COMMENT" &> "KNN_$FILE_SUFFIX.out" &
python -u -m cProfile -o "NB_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt NaiveBayes -d "NB $COMMENT" &> "NB_$FILE_SUFFIX.out" &
# python -u -m cProfile -o "BLR_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt BayesianLogisticRegression -d "BLR $COMMENT" &> "BLR_$FILE_SUFFIX.out" &