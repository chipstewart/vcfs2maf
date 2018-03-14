#!/bin/bash 
echo $@
PID=$1
MAF1=$2  
NCALL=$3
ALGS=$4
export JULIA_LOAD_PATH="/opt/src"
echo "julia /opt/src/rebc_consensus_maf.jl $MAF1 $NCALL $ALGS $PID"
julia /opt/src/rebc_consensus_maf.jl $MAF1  $NCALL $ALGS $PID


