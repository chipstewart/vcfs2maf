VCF1=$3  # M1 format
VCF2=$4  # M2 format vcf
VCF3=$5  # Strelka SNV format
VCF4=$6  # Strelka indel format
TID=$1
NID=$2
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/vcf2txt.jl $VCF1  tmp1.tsv  
julia /opt/src/vcf2txt.jl $VCF2  tmp2.tsv
julia /opt/src/vcf2txt.jl $VCF3  tmp3.tsv
julia /opt/src/vcf2txt.jl $VCF4  tmp4.tsv
julia /opt/src/StrelkaM1M2_maflite.jl $TID $NID tmp1.tsv tmp2.tsv tmp3.tsv tmp4.tsv m1m2strelka_maflite.tsv


