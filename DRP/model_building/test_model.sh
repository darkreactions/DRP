#!/usr/bin/env bash
source DRP/model_building/set_class_path.sh
TMP_DIR="tmp"

java weka.classifiers.trees.J48 -T $TMP_DIR/$1_out.arff -l $2 -p 0 -c last > $TMP_DIR/$1.out >&2

#java weka.classifiers.trees.J48 -T $1 -l $2 -p 0 -c last 1> $3

