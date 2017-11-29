PID=$1
INMAF=$2
REMOVE=$3
OUTMAF=$4  # Strelka SNV format
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/trim_maf_fields.jl $PID $INMAF  $REMOVE $PID.annotated.trim.maf


