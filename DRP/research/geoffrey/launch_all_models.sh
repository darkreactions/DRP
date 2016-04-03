#!/bin/bash
SCRIPT="../build_model.py"
DESCRIPTOR_FILE="descs/new_legacy_newCA_nonZeroVariance.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
COMMENT="BCR Weighted. Legacy rxns. New and Legacy descriptors new ChemAxon only NonzeroVariance. 15 mutual info split."
FILE_SUFFIX="BCR_legRxnNonZeroCompound_NewLegDscNewCANonzeroVar_15MISplit"
SPLITTER="MutualInfoSplitter"

#python -u -m cProfile -o "SVMbasic_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_basic -d "SVM Default PUK $COMMENT" &> "SVMbasic_$FILE_SUFFIX.out" &
python -u -m cProfile -o "J48_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt J48 -d "J48 $COMMENT" &> "J48_$FILE_SUFFIX.out" &
python -u -m cProfile -o "KNN_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt KNN -d "KNN $COMMENT" &> "KNN_$FILE_SUFFIX.out" &
python -u -m cProfile -o "NB_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt NaiveBayes -d "NB $COMMENT" &> "NB_$FILE_SUFFIX.out" &
python -u -m cProfile -o "RF_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt RandomForest -d "Random Forest $COMMENT" &> "RF_$FILE_SUFFIX.out" &
python -u -m cProfile -o "LR_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt LogisticRegression -d "Logistic Regression $COMMENT" &> "LR_$FILE_SUFFIX.out" &
# python -u -m cProfile -o "BLR_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt BayesianLogisticRegression -d "BLR $COMMENT" &> "BLR_$FILE_SUFFIX.out" &