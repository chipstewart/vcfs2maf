PID=$1
TID=$2
NID=$3
TBAM=$4
NBAM=$5
VCF_INDEL=$6  # SvABA indel format
MAX_NORMAL_ALT_COUNT=$7
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/vcf2txt.jl $VCF_INDEL  $PID.SvABA.INDEL.tsv
julia /opt/src/svaba_maflite.jl $TID $NID  $TBAM $NBAM  $PID.SvABA.INDEL.tsv $PID.SvABA_maflite.tsv $MAX_NORMAL_ALT_COUNT


