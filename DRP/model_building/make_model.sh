#!/usr/bin/env bash
source DRP/model_building/set_class_path.sh

<<<<<<< HEAD
# PUK Kernel:
#java weka.classifiers.functions.SMO -d $1 -t $2 -K "weka.classifiers.functions.supportVector.Puk -O 0.5 -S 7" -p 0

# Linear Kernel (Default)
java weka.classifiers.functions.SMO -d $1 -t $2 -p 0
=======
java weka.classifiers.functions.SMO -d $1 -t $2 -K "weka.classifiers.functions.supportVector.Puk -O 0.5 -S 7" -p 0   # PUK Kernel

#java weka.classifiers.functions.SMO -d $1 -t $2 -p 0   # Linear Kernel
>>>>>>> e914b32092c53d44b4feb0084f9db053af498f4d
