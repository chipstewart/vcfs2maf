PID=$1
TID=$2
NID=$3
VCF_SNV=$4  # Strelka SNV format
VCF_INDEL=$5  # Strelka indel format
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/vcf2txt.jl $VCF_SNV  $PID.StrelkaSNV.tsv  
julia /opt/src/vcf2txt.jl $VCF_INDEL  $PID.StrelkaINDEL.tsv
julia /opt/src/Strelka_maflite.jl $TID $NID  $PID.StrelkaSNV.tsv $PID.StrelkaINDEL.tsv $PID.Strelka_maflite.tsv


