VCF1=$3
VCF2=$4
TID=$1
NID=$2
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/vcf2txt.jl $VCF1  tmp1.tsv
julia /opt/src/vcf2txt.jl $VCF2  tmp2.tsv
julia /opt/src/M1M2maflite.jl $TID $NID tmp1.tsv tmp2.tsv m1m2_maflite.tsv


