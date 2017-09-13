TID=$1
NID=$2
PID=$3
VCF1=$4
export JULIA_LOAD_PATH="/opt/src"

julia /opt/src/vcf2txt.jl $VCF1  tmp1.tsv
sed "s/${TID}/TUMOR/g" tmp1.tsv > tmp2.tsv
sed "s/${NID}/NORMAL/g" tmp2.tsv > tmp3.tsv


julia /opt/src/gatk4_m2_maflite.jl $TID $NID tmp3.tsv $PID


