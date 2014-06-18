#!/usr/bin/env bash
source set_class_path.sh

java weka.classifiers.trees.J48 -T $1 -l $2 -p 0 -c last > $3

