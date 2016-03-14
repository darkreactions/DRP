SCRIPT="../build_model.py"
DESCRIPTOR_FILE="descs/legacy_valid_allNewReplacements_removeGroupPeriodValenceSlowCoolLeak.dsc"
OUTCOME_DESCRIPTOR="boolean_outcome_legacy"
RXN_SET_NAME="valid_legacy_reactions"
COMMENT="Legacy descriptors, all new direct replacements - Group, Period, Valence, Slow Cool, Leak. 15 mutual info splits."
FILE_SUFFIX="legRxnNonZeroCompound_legDscValidAllNewReplacementsRemoveGroupPeriodValenceSlowCoolLeak"
SPLITTER="MutualInfoSplitter"

python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_basic -d "SVM $COMMENT" &> "SVMbasic_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt SVM_PUK_BCR -d "$BCR SVM $COMMENT" &> "SVMBCR_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt J48 -d "J48 $COMMENT" &> "J48_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt KNN -d "KNN $COMMENT" &> "KNN_$FILE_SUFFIX.out" &
python -u $SCRIPT -p @$DESCRIPTOR_FILE -rxn $RXN_SET_NAME -r $OUTCOME_DESCRIPTOR -s $SPLITTER -v -mt NaiveBayes -d "NB $COMMENT" &> "NB_$FILE_SUFFIX.out" &
