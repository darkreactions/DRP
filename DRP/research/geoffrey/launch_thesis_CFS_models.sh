python -u build_model.py -p @descs/CFS_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt SVM_PUK_basic -d "geoffrey thesis SVM CFS" &> SVM_basic_CFS_thesis.out &
python -u build_model.py -p @descs/CFS_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt SVM_PUK_BCR -d "geoffrey thesis BCR SVM CFS" &> SVM_BCR_CFS_thesis.out &
python -u build_model.py -p @descs/CFS_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt KNN -d "geoffrey thesis BCR SVM CFS" &> KNN_CFS_thesis.out &
python -u build_model.py -p @descs/CFS_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt J48 -d "geoffrey thesis BCR SVM CFS" &> J48_CFS_thesis.out &
python -u build_model.py -p @descs/CFS_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt NaiveBayes -d "geoffrey thesis BCR SVM CFS" &> NB_CFS_thesis.out &
