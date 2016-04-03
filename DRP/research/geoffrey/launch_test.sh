#!/bin/bash
SCRIPT="../build_many_models.py"
DESCRIPTOR_FILE="descs/new_noCA_nonZeroInfo.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_rxns_nonzero_compound"
COMMENT="TEST"
FILE_SUFFIX="legRxnNonZeroCompound_NewNoCANonzeroInfo_15MISplit"
SPLITTER="MutualInfoSplitter"
SPLITTER_OPTIONS="{'num_splits': 15}"
MODEL_VISITOR_FILE="mv.txt"
MODEL_VISITOR_OPTION_FILE="vo.txt"

python $SCRIPT -p 'test' -s $SPLITTER