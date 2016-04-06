#!/bin/bash

SCRIPT="../build_many_models.py"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
SPLITTER="RandomSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
MODEL_VISITOR_FILE="mv.txt"
MODEL_VISITOR_OPTION_FILE="vo.txt"

DESCRIPTOR_DIR="descs"
DESCRIPTOR_FILES=(
    #"legacy_nonZeroVariance.dsc"
    #"legacy_noCA_nonZeroVariance.dsc"
    #"new_noCA_nonZeroVariance.dsc"
    #"new_CA_nonZeroVariance.dsc"
    #"new_legacy_CA_nonZeroVariance.dsc"
    #"new_legacy_newCA_nonZeroVariance.dsc"
    #"new_legacy_noCA_nonZeroVariance.dsc"
    #"new_legacy_legCA_nonZeroVariance.dsc"
    #"new_reactionpH_nonZeroVariance.dsc"

    "legacy_noCA_nonZeroInfo.dsc"
    "legacy_nonZeroInfo.dsc"
    "new_CA_nonZeroInfo.dsc"
    "new_legacy_CA_nonZeroInfo.dsc"
    # new_legacy_legCA_nonZeroInfo.dsc
    # new_legacy_newCA_nonZeroInfo.dsc
    # new_legacy_noCA_nonZeroInfo.dsc
    # new_noCA_nonZeroInfo.dsc
    # new_reactionpH_nonZeroInfo.dsc
)

for DESCRIPTOR_FN in "${DESCRIPTOR_FILES[@]}"
do
    DESCRIPTOR_FILE="$DESCRIPTOR_DIR/$DESCRIPTOR_FN"
    COMMENT="Legacy rxns. ${DESCRIPTOR_FN}. 15 Random split."
    FILE_SUFFIX="models_${DESCRIPTOR_FN}_15RandSplit"

    python -u $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -v -mt @"$MODEL_VISITOR_FILE" -vo @"$MODEL_VISITOR_OPTION_FILE" -d "$COMMENT" &> "$FILE_SUFFIX.out" &
done