#!/bin/bash
SCRIPT="../run_feature_selection.py"
# legacy_noCA_nonZeroVariance.dsc
# new_CA_nonZeroVariance.dsc
# new_legacy_legCA_nonZeroVariance.dsc
# new_legacy_noCA_nonZeroVariance.dsc

DESCRIPTOR_FILE="descs/new_CA_nonZeroVariance.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
COMMENT="Legacy rxns. New descriptors with ChemAxon only NonzeroVariance."
FILE_SUFFIX="legRxnNonZeroCompound_NewDscCANonzeroVar"
SPLITTER="MutualInfoSplitter"

python -u -m cProfile -o "CFS_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -trs $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -v -vt CFS -d "CFS $COMMENT" -o "CFS_$FILE_SUFFIX.dsc" &> "CFS_$FILE_SUFFIX.out" &
#python -u -m cProfile -o "InfoGain_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -trs $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -v -vt InfoGain -d "InfoGain $COMMENT" -o "InfoGain_$FILE_SUFFIX.dsc" &> "InfoGain_$FILE_SUFFIX.out" &
#python -u -m cProfile -o "ChiSquared_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -trs $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -v -vt ChiSquared -d "ChiSquared $COMMENT" -o "ChiSquared_$FILE_SUFFIX.dsc" &> "ChiSquared_$FILE_SUFFIX.out" &
