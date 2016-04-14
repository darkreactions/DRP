#!/bin/bash

SCRIPT="../build_many_models.py"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
SPLITTER="MutualInfoSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
MODEL_VISITOR_FILE="mv_restart.txt"
MODEL_VISITOR_OPTION_FILE="vo_restart.txt"


DESCRIPTOR_DIR="final_descs/use"
DESCRIPTOR_FILES=(
    # Run:
    # "new_legacy_newCA_nonZeroInfo.dsc"
    # "new_leak_slowcool_group_period_valence_nonZeroInfo.dsc"
    # "new_legacy_legCA_noPSA_nonZeroVariance.dsc"
    # "new_legacy_legCA_noPSA_nonZeroInfo.dsc"
    # "new_legacy_bothCA_noPSA_nonZeroInfo.dsc"
    # "new_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
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

    # Running:
    "new_legacy_newCA_nonZeroVariance.dsc"


    # Not yet Run:


)

for DESCRIPTOR_FN in "${DESCRIPTOR_FILES[@]}"
do
    DESCRIPTOR_FILE="$DESCRIPTOR_DIR/$DESCRIPTOR_FN"
    COMMENT="Legacy rxns. ${DESCRIPTOR_FN}. 15 MI split."
    FILE_SUFFIX="models_${DESCRIPTOR_FN}_15MISplit"

    python -u $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -v -mt @"$MODEL_VISITOR_FILE" -vo @"$MODEL_VISITOR_OPTION_FILE" -d "$COMMENT" &>> "$FILE_SUFFIX.out" &
done