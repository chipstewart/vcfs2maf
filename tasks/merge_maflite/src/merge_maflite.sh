PID=$1
TID=$2
NID=$3
MAF1=$4  # M1+M2
MAF2=$5  # Strelka
MAF3=$6  # SvABA 
MAF4=$7  # Other
LAB1=$8  # M1M2
LAB2=$9  # STRELKA
LAB3=$10  # SVABA
LAB4=$11  # OTHER
export JULIA_LOAD_PATH="/opt/src"
julia /opt/src/merge_maflite.jl $TID $NID $MAF1 $MAF2 $MAF3 $MAF4 $LAB1 $LAB2  $LAB3  $LAB4 $PID.merged.maflite.tsv


