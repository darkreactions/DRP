#!/bin/bash
SCRIPT="../run_feature_selection.py"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
FILE_SUFFIX="legRxnNonZeroCompound_NewDscNoCANonzeroVar"

DESCRIPTOR_DIR="final_descs/use"
DESCRIPTOR_FILES=(
    "legacy_mw_noCA_nonZeroVariance.dsc"
    "new_legacy_noCA_nonZeroVariance.dsc"
    "legacy_mw_noPSA_nonZeroVariance.dsc"
    "new_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
    "new_legacy_bothCA_noPSA_nonZeroVariance.dsc"
    "new_legacy_legCA_noPSA_nonZeroVariance.dsc"
    "new_legacy_newCA_nonZeroVariance.dsc"
    "new_noCA_leak_slowcool_group_period_valence_nonZeroVariance.dsc"
)

for DESCRIPTOR_FN in "${DESCRIPTOR_FILES[@]}"
do
    DESCRIPTOR_FILE="$DESCRIPTOR_DIR/$DESCRIPTOR_FN"
    COMMENT="Legacy rxns. ${DESCRIPTOR_FN}"
    FILE_SUFFIX="${DESCRIPTOR_FN}"

    python -u $SCRIPT -p @$DESCRIPTOR_FILE -trs $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -v -vt CFS -d "CFS $COMMENT" -o "CFS_$FILE_SUFFIX.dsc" &> "CFS_$FILE_SUFFIX.out"
done



#python -u -m cProfile -o "InfoGain_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -trs $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -v -vt InfoGain -d "InfoGain $COMMENT" -o "InfoGain_$FILE_SUFFIX.dsc" &> "InfoGain_$FILE_SUFFIX.out" &
#python -u -m cProfile -o "ChiSquared_$FILE_SUFFIX.profile" $SCRIPT -p @$DESCRIPTOR_FILE -trs $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -v -vt ChiSquared -d "ChiSquared $COMMENT" -o "ChiSquared_$FILE_SUFFIX.dsc" &> "ChiSquared_$FILE_SUFFIX.out" &
