#!/usr/bin/env bash
source set_class_path.sh
$TMP_DIR="/home/drp/web/darkreactions.haverford.edu/app/DRP/tmp"


java weka.classifiers.trees.J48 -T $TMP_DIR/$1_out.arff -l /home/drp/research/2.3.14_svm_true.model -p 0 -c last > $TMP_DIR/$1.out


