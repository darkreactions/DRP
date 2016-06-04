#!/bin/bash

SCRIPT="../split_dump_arffs.py"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
SPLITTER="MutualInfoSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"

DESCRIPTOR_DIR="../legacy_tests/final_descs/use"
DESCRIPTOR_FILES=(
    # Run:

    # Running:

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
    # "new_legacy_newCA_nonZeroVariance.dsc"

    "CFS_new_legacy_bothCA_noPSA_nonZeroVariance.dsc"
    "CFS_new_legacy_legCA_noPSA_nonZeroVariance.dsc"
    "CFS_new_legacy_newCA_nonZeroVariance.dsc"
    "CFS_new_legacy_noCA_nonZeroVariance.dsc"



    # Not yet Run:


)

for DESCRIPTOR_FN in "${DESCRIPTOR_FILES[@]}"
do
    DESCRIPTOR_FILE="$DESCRIPTOR_DIR/$DESCRIPTOR_FN"
    COMMENT="${DESCRIPTOR_FN}.15_MISplit."
    FILE_SUFFIX="splitting_15MI_${DESCRIPTOR_FN}"

    python -u $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -v -d "$COMMENT" &>> "$FILE_SUFFIX.out"
done