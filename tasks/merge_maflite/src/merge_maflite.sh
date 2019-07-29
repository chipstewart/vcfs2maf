#!/bin/bash 
echo $@
PID=$1
TID=$2
NID=$3
MAF1=$4  # M1
MAF2=$5  # M2
MAF3=$6  # Strelka 
MAF4=$7  # Snowman
MAF5=$8  # SvABA 
MAF6=$9  # Strelka2
LAB1=${10}  # M1
LAB2=${11}  # M2
LAB3=${12}  # STRELKA
LAB4=${13}  # Snowman
LAB5=${14}  # SVABA
LAB6=${15}  # Strelka2
BUILD=${16} # reference build 
export JULIA_LOAD_PATH="/opt/src"
echo "julia /opt/src/merge_maflite.jl $TID $NID $MAF1 $MAF2 $MAF3 $MAF4 $MAF5 $MAF6 $LAB1 $LAB2 $LAB3 $LAB4 $LAB5 $LAB6 $PID.merged.maflite.tsv $BUILD" 
julia /opt/src/merge_maflite.jl $TID $NID $MAF1 $MAF2 $MAF3 $MAF4 $MAF5 $MAF6 $LAB1 $LAB2 $LAB3 $LAB4 $LAB5 $LAB6 $PID.merged.maflite.tsv $BUILD


