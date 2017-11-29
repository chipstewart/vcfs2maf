#!/bin/bash 
echo $@
PID=$1
TID=$2
NID=$3
MAF1=$4  # M1
MAF2=$5  # M2
MAF3=$6  # Strelka 
MAF4=$7  # SvABA 
MAF5=$8  # SvABA 
LAB1=$9  # M1
LAB2=${10}  # M2
LAB3=${11}  # STRELKA
LAB4=${12}  # SVABA
LAB5=${13}  # SVABA
export JULIA_LOAD_PATH="/opt/src"
echo "julia /opt/src/merge_maflite.jl $TID $NID $MAF1 $MAF2 $MAF3 $MAF4 $MAF5 $LAB1 $LAB2  $LAB3  $LAB4 $PID.merged.maflite.tsv"
julia /opt/src/merge_maflite.jl $TID $NID $MAF1 $MAF2 $MAF3 $MAF4 $MAF5 $LAB1 $LAB2  $LAB3  $LAB4  $LAB5 $PID.merged.maflite.tsv


