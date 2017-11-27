TID=$1
NID=$2
PID=$3
MAFLITE1=$4
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/m1_maflite.jl $TID $NID $MAFLITE1 $PID.m1_maflite.tsv


