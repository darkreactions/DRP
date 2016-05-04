#!/bin/bash

SCRIPT="../validate_model.py"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="new_rxns_nonzero_compound"
SPLITTER="RandomSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
MODEL_TOOL="J48"
MODEL_VISITOR_OPTIONS="{'BCR':False}"

DESCRIPTOR_DIR="final_descs/use"
DESCRIPTOR_FILES=(
    # "new_legacy_newCA_nonZeroInfo.dsc"
    # "new_leak_slowcool_group_period_valence_nonZeroInfo.dsc"
    # "new_legacy_legCA_noPSA_nonZeroVariance.dsc"
    # "new_legacy_legCA_noPSA_nonZeroInfo.dsc"
    # "new_legacy_bothCA_noPSA_nonZeroInfo.dsc"
    "new_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
    # "legacy_mw_noPSA_nonZeroInfo.dsc"
    # "legacy_mw_noPSA_nonZeroVariance.dsc"
    # "new_legacy_noCA_nonZeroVariance.dsc"
    # "new_noCA_leak_slowcool_group_period_valence_nonZeroInfo.dsc"
    # "new_noCA_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
    # "new_legacy_noCA_nonZeroInfo.dsc"
    # "new_legacy_bothCA_noPSA_nonZeroVariance.dsc"
    # "legacy_mw_noCA_nonZeroVariance.dsc"
    # "legacy_mw_noCA_nonZeroInfo.dsc"
    # "CFS_legacy_mw_noCA_nonZeroVariance.dsc"
    # "CFS_legacy_mw_noPSA_nonZeroVariance.dsc"
    # "CFS_new_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
    # "CFS_new_noCA_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
    # "new_legacy_newCA_nonZeroVariance.dsc"
)

for DESCRIPTOR_FN in "${DESCRIPTOR_FILES[@]}"
do
    DESCRIPTOR_FILE="$DESCRIPTOR_DIR/$DESCRIPTOR_FN"
    FILE_SUFFIX="predictions_from_${MODEL_TOOL}_${DESCRIPTOR_FN}_15RandSplit"

    python -u $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -v -mt "$MODEL_TOOL" -vo "$MODEL_VISITOR_OPTIONS" &>> "$FILE_SUFFIX.out" &
done