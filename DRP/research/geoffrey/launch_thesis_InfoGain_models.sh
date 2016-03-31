python -u build_model.py -p @descs/InfoGain_top10_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt SVM_PUK_basic -d "geoffrey thesis SVM InfoGain_top10" &> SVM_basic_InfoGain_top10_thesis.out &
python -u build_model.py -p @descs/InfoGain_top10_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt SVM_PUK_BCR -d "geoffrey thesis BCR SVM InfoGain_top10" &> SVM_BCR_InfoGain_top10_thesis.out &
python -u build_model.py -p @descs/InfoGain_top10_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt KNN -d "geoffrey thesis BCR SVM InfoGain_top10" &> KNN_InfoGain_top10_thesis.out &
python -u build_model.py -p @descs/InfoGain_top10_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt J48 -d "geoffrey thesis BCR SVM InfoGain_top10" &> J48_InfoGain_top10_thesis.out &
python -u build_model.py -p @descs/InfoGain_top10_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt NaiveBayes -d "geoffrey thesis BCR SVM InfoGain_top10" &> NB_InfoGain_top10_thesis.out &
