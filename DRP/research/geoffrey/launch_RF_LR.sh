#!/bin/bash

SCRIPT="../build_many_models.py"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
SPLITTER="MutualInfoSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
MODEL_VISITOR_FILE="mv_RF_LR.txt"
MODEL_VISITOR_OPTION_FILE="vo_RF_LR.txt"

DESCRIPTOR_DIR="descs"
DESCRIPTOR_FILES=(
    # "legacy_noCA_nonZeroVariance.dsc"
    # "new_noCA_nonZeroVariance.dsc"
    # "new_CA_nonZeroVariance.dsc"
    "new_legacy_CA_nonZeroVariance.dsc"
    "new_legacy_newCA_nonZeroVariance.dsc"
    "new_legacy_noCA_nonZeroVariance.dsc"
    "new_reactionpH_nonZeroVariance.dsc"
)

for DESCRIPTOR_FN in "${DESCRIPTOR_FILES[@]}"
do
    DESCRIPTOR_FILE="$DESCRIPTOR_DIR/$DESCRIPTOR_FN"
    COMMENT="Legacy rxns. ${DESCRIPTOR_FN}. 15 MI split."
    FILE_SUFFIX="RF_LR_${DESCRIPTOR_FN}_15MISplit"

    python -u $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -v -mt @"$MODEL_VISITOR_FILE" -vo @"$MODEL_VISITOR_OPTION_FILE" -d "$COMMENT" &> "$FILE_SUFFIX.out" &
done