#!/usr/bin/env bash
source DRP/model_building/set_class_path.sh

# PUK Kernel:
#java weka.classifiers.functions.SMO -d $1 -t $2 -K "weka.classifiers.functions.supportVector.Puk -O 0.5 -S 7" -p 0

# Linear Kernel (Default)
java weka.classifiers.functions.SMO -d $1 -t $2 -p 0
