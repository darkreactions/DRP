#!/usr/bin/env bash
source DRP/model_building/set_class_path.sh
TMP_DIR="tmp"

java weka.classifiers.trees.J48 -T $TMP_DIR/$1_out.arff -l /home/drp/research/2.3.14_svm_true.model -p 0 -c last > $TMP_DIR/$1.out >&2


