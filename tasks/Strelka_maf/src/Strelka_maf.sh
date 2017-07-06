VCF_SNV=$1  # Strelka SNV format
VCF_INDEL=$2  # Strelka indel format
TID=$1
NID=$2
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/vcf2txt.jl $VCF_SNV  StrelkaSNV.tsv  
julia /opt/src/vcf2txt.jl $VCF_INDEL  StrelkaINDEL.tsv
julia /opt/src/Strelka_maflite.jl $TID $NID  StrelkaSNV.tsv StrelkaINDEL.tsv Strelka_maflite.tsv


