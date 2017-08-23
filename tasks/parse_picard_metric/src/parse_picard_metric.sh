#!/bin/bash 
echo $@
FILE=$1
FIELD=$2  
export JULIA_LOAD_PATH="/opt/src"
echo "julia /opt/src/parse_picard_metric.jl $FILE $FIELD "
julia /opt/src/parse_picard_metric.jl $FILE $FIELD
