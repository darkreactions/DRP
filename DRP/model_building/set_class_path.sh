#!/usr/bin/env bash

# Get the weka path specified by data_config.
WEKA_PATH=`python DRP/data_config.py weka_path`

# Add that path to the classpath.
export CLASSPATH=$CLASSPATH:$WEKA_PATH

