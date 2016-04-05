#!/bin/bash
SCRIPT="../build_many_models.py"
# 

DESCRIPTOR_FILE="descs/CFS_legRxnNonZeroCompound_NewLegDscNoCANonzeroVar.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
COMMENT="Legacy rxns. CFS of new and legacy descriptors without ChemAxon. 15 mutual info split."
FILE_SUFFIX="legRxnNonZeroCompound_CFSNewLegDscNoCA_15MISplit"
SPLITTER="MutualInfoSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
MODEL_VISITOR_FILE="mv.txt"
MODEL_VISITOR_OPTION_FILE="vo.txt"

python -u -m cProfile -o "$FILE_SUFFIX.profile" $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -v -mt @"$MODEL_VISITOR_FILE" -vo @"$MODEL_VISITOR_OPTION_FILE" -d "$COMMENT" &> "$FILE_SUFFIX.out" &
