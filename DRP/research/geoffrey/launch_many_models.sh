#!/bin/bash

DESCRIPTOR_FILE="descs/CFS_legRxnNonZeroCompound_NewLegDscLegCANonzeroVar.dsc"
COMMENT="Legacy rxns. CFS of new and legacy descriptors with legacy ChemAxon. 15 mutual info split."
FILE_SUFFIX="restart_legRxnNonZeroCompound_CFSNewLegDscLegCA_15MISplit"


SCRIPT="../build_many_models.py"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
SPLITTER="MutualInfoSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
MODEL_VISITOR_FILE="mv_restart.txt"
MODEL_VISITOR_OPTION_FILE="vo_restart.txt"

python -u -m cProfile -o "$FILE_SUFFIX.profile" $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -v -mt @"$MODEL_VISITOR_FILE" -vo @"$MODEL_VISITOR_OPTION_FILE" -d "$COMMENT" &> "$FILE_SUFFIX.out" &
