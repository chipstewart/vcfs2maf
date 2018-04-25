PID=$1
INMAF=$2
COPY_FIELDS=$3
EXTENSION=$4  # Strelka SNV format
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/copy_maf_fields.jl $PID $INMAF  $COPY_FIELDS $EXTENSION


