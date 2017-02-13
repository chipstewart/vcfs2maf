VCF1=$1
VCF2=$2
TID=$3
NID=$4
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/vcf2txt.jl $VCF1  ./tmp1.tsv
julia /opt/src/vcf2txt.jl $VCF2  ./tmp2.tsv
julia /opt/src/M1M2maflite.jl $TID $NID ./test1.tsv  ./test2.tsv m1m2_maflite.tsv
tar cvfz tmp1.tsv tmp1.tsv 


