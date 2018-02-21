#!/bin/bash 
echo $@
PID=$1
MAF1=$2  
CONTEST=$3
PTHRESH=$4
export JULIA_LOAD_PATH="/opt/src"
echo "julia /opt/src/postFilter_contEst.jl $PID $MAF1 $CONTEST $PTHRESH "
julia /opt/src/postFilter_contEst.jl  $PID $MAF1 $CONTEST $PTHRESH 


