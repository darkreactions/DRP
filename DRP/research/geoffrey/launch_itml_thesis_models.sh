python -u build_model.py -p @descs/itml_10k_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt SVM_PUK_basic -d "geoffrey thesis SVM itml 10k descriptors, no feature selection" &> SVM_basic_itml_10k_thesis.out &
python -u build_model.py -p @descs/itml_10k_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt SVM_PUK_BCR -d "geoffrey thesis BCR SVM itml 10k descriptors, no feature selection" &> SVM_BCR_itml_10k_thesis.out &
python -u build_model.py -p @descs/itml_10k_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt KNN -d "geoffrey thesis BCR SVM itml 10k descriptors, no feature selection" &> KNN_itml_10k_thesis.out &
python -u build_model.py -p @descs/itml_10k_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt J48 -d "geoffrey thesis BCR SVM itml 10k descriptors, no feature selection" &> J48_itml_10k_thesis.out &
python -u build_model.py -p @descs/itml_10k_thesis.dsc -trs geoffrey_split_for_thesis_0 -tes geoffrey_split_for_thesis_1 -v -mt NaiveBayes -d "geoffrey thesis BCR SVM itml 10k descriptors, no feature selection" &> NB_itml_10k_thesis.out &
