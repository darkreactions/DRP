#!/bin/bash

SCRIPT="../optimize_PUK.py"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
SPLITTER="MutualInfoSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
TYPE=$1

SIGMA_MIN=0.25
SIGMA_MAX=64
SIGMA_STEP=2
OMEGA_MIN=0.25
OMEGA_MAX=64
OMEGA_STEP=2



DESCRIPTOR_DIR="final_descs/use"
DESCRIPTOR_FILES=(
    # "CFS_legacy_mw_noCA_nonZeroVariance.dsc"
    # "CFS_legacy_mw_noPSA_nonZeroVariance.dsc"
    # "CFS_new_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
    # "CFS_new_legacy_bothCA_noPSA_nonZeroVariance.dsc"
    # "CFS_new_legacy_legCA_noPSA_nonZeroVariance.dsc"
    # "CFS_new_legacy_newCA_nonZeroVariance.dsc"
    # "CFS_new_legacy_noCA_nonZeroVariance.dsc"
    # "CFS_new_noCA_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
    "legacy_mw_noCA_nonZeroInfo.dsc"
    # "legacy_mw_noCA_nonZeroVariance.dsc"
    # "legacy_mw_noPSA_nonZeroInfo.dsc"
    # "legacy_mw_noPSA_nonZeroVariance.dsc"
    # "new_leak_slowcool_group_period_valence_nonZeroInfo.dsc"
    # "new_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
    # "new_legacy_bothCA_noPSA_nonZeroInfo.dsc"
    # "new_legacy_bothCA_noPSA_nonZeroVariance.dsc"
    # "new_legacy_legCA_noPSA_nonZeroInfo.dsc"
    # "new_legacy_legCA_noPSA_nonZeroVariance.dsc"
    # "new_legacy_newCA_nonZeroInfo.dsc"
    # "new_legacy_newCA_nonZeroVariance.dsc"
    "new_legacy_noCA_nonZeroInfo.dsc"
    # "new_legacy_noCA_nonZeroVariance.dsc"
    "new_noCA_leak_slowcool_group_period_valence_nonZeroInfo.dsc"
    # "new_noCA_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
)



if [ "$TYPE" = "W" ] ; then
    VISITOR_OPTIONS="{'BCR': True}"
else
    if [ "$TYPE" = "U" ] ; then
	VISITOR_OPTIONS="{}"
    else
	echo "TYPE must be weighted or unweighted"
	exit
    fi
fi

for DESCRIPTOR_FN in "${DESCRIPTOR_FILES[@]}"
do
    DESCRIPTOR_FILE="$DESCRIPTOR_DIR/$DESCRIPTOR_FN"
    COMMENT="Legacy rxns. Optimizing PUK. $TYPE ${DESCRIPTOR_FN}. 15 MI split."
    FILE_SUFFIX="PUK_scan_${TYPE}_${DESCRIPTOR_FN}_15MISplit"

    python -u $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -vo "$VISITOR_OPTIONS" -v -d "$COMMENT" -ps $SIGMA_MIN $SIGMA_MAX $SIGMA_STEP -po $OMEGA_MIN $OMEGA_MAX $OMEGA_STEP -g &> "$FILE_SUFFIX.out" &
done
