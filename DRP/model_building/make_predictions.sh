#!/usr/bin/env bash
source DRP/model_building/set_class_path.sh

java weka.classifiers.trees.J48 -T $1 -l $2 -p 0 -c last 1> $3

