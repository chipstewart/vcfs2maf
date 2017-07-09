TID=$1
NID=$2
PID=$3
VCF1=$4
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/vcf2txt.jl $VCF1  tmp1.tsv
julia /opt/src/m2_maflite.jl $TID $NID tmp1.tsv $PID.m2_maflite.tsv


