PID=$1
TID=$2
NID=$3
VCF_SNV=$4  # Strelka SNV format
VCF_INDEL=$5  # Strelka indel format
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/vcf2txt.jl $VCF_SNV  $PID.Strelka2SNV.tsv  
julia /opt/src/vcf2txt.jl $VCF_INDEL  $PID.Strelka2INDEL.tsv
julia /opt/src/strelka2_maflite.jl $TID $NID  $PID.Strelka2SNV.tsv $PID.Strelka2INDEL.tsv $PID.Strelka2_maflite.tsv


