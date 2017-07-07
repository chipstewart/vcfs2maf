M1M2MAF=$3
SMVCF=$4
TID=$1
NID=$2
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/vcf2txt.jl $SMVCF  SM.tsv
julia /opt/src/M1M2mafliteSnowman.jl $TID $NID $M1M2MAF SM.tsv m1m2SM_maflite.tsv


