PID=$1
TID=$2
NID=$3
MAF1=$4  # M1
MAF2=$5  # M2
MAF3=$6  # Strelka 
MAF4=$7  # SvABA 
LAB1=$8  # M1
LAB2=$9  # M2
LAB3=$10  # STRELKA
LAB4=$11  # SVABA
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/merge_maflite.jl $TID $NID $MAF1 $MAF2 $MAF3 $MAF4 $LAB1 $LAB2  $LAB3  $LAB4 $PID.merged.maflite.tsv


