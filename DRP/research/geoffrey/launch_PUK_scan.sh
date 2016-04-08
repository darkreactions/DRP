#!/bin/bash

SCRIPT="../optimize_PUK.py"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
SPLITTER="MutualInfoSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
TYPE=$1

SIGMA_MIN=0.25
SIGMA_MAX=64
SIGMA_STEP=4
OMEGA_MIN=0.25
OMEGA_MAX=64
OMEGA_STEP=4



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

    # "legacy_noCA_nonZeroInfo.dsc"
    # "legacy_nonZeroInfo.dsc"
    # "new_reactionpH_nonZeroInfo.dsc"
    # "new_CA_nonZeroInfo.dsc"
    # "new_legacy_noCA_nonZeroInfo.dsc"
    # "new_legacy_CA_nonZeroInfo.dsc"
    # "new_legacy_legCA_nonZeroInfo.dsc"
    # "new_legacy_newCA_nonZeroInfo.dsc"
    # "new_noCA_nonZeroInfo.dsc"

    "CFS_legRxnNonZeroCompound_LegDscCANonzeroVar.dsc"
    # "CFS_legRxnNonZeroCompound_LegDscNoCANonzeroVar.dsc"
    # "CFS_legRxnNonZeroCompound_NewDscCANonzeroVar.dsc"
    # "CFS_legRxnNonZeroCompound_NewDscNoCANonzeroVar.dsc"
    # "CFS_legRxnNonZeroCompound_NewDscRxnpHCANonzeroVar.dsc"
    # "CFS_legRxnNonZeroCompound_NewLegDscBothCANonzeroVar.dsc"
    # "CFS_legRxnNonZeroCompound_NewLegDscLegCANonzeroVar.dsc"
    # "CFS_legRxnNonZeroCompound_NewLegDscNewCANonzeroVar.dsc"
    # "CFS_legRxnNonZeroCompound_NewLegDscNoCANonzeroVar.dsc"
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
    FILE_SUFFIX="PUK_SCAN_${TYPE}_${DESCRIPTOR_FN}_15MISplit"

    python -u $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -vo "$VISITOR_OPTIONS" -v -d "$COMMENT" -ps $SIGMA_MIN $SIGMA_MAX $SIGMA_STEP -po $OMEGA_MIN $OMEGA_MAX $OMEGA_STEP -g &> "$FILE_SUFFIX.out" &
done
