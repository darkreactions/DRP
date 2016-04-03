#!/bin/bash
SCRIPT="../build_many_models.py"
DESCRIPTOR_FILE="descs/legacy_noCA_nonZeroInfo.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
COMMENT="Legacy rxns. Legacy descriptors without ChemAxon nonzero InfoGain. 15 mutual info split."
FILE_SUFFIX="legRxnNonZeroCompound_LegNoCANonzeroInfo_15MISplit"
SPLITTER="MutualInfoSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
MODEL_VISITOR_FILE="mv.txt"
MODEL_VISITOR_OPTION_FILE="vo.txt"

python -u -m cProfile -o "$FILE_SUFFIX.profile" $SCRIPT -p @"$DESCRIPTOR_FILE" -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -so "$SPLITTER_OPTIONS" -v -mt @"$MODEL_VISITOR_FILE" -vo @"$MODEL_VISITOR_OPTION_FILE" -d "$COMMENT" &> "$FILE_SUFFIX.out" &
