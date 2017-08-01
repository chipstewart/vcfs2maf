#!/bin/bash 
echo $@
PID=$1
MAF1=$2  
NCALL=$3
LONGINDEL=$4
export JULIA_LOAD_PATH="/opt/src"
echo "julia /opt/src/consensus_maf.jl $MAF1 $NCALL $LONGINDEL $PID.consensus.maf"
julia /opt/src/consensus_maf.jl $MAF1  $NCALL $LONGINDEL $PID.consensus.maf


