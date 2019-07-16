TID=$1
NID=$2
PID=$3
VCF1=$4
BUILD=$5
DETIN_FILT=$6

set -x

echo $TID
echo $NID
echo $PID
echo $VCF1
echo $BUILD
echo $DETIN_FILT

export JULIA_LOAD_PATH=":.:/opt/src"
which julia

julia /opt/src/vcf2txt.jl $VCF1  tmp1.tsv
sed "s/${TID}/TUMOR/g" tmp1.tsv > tmp2.tsv
sed "s/${NID}/NORMAL/g" tmp2.tsv > tmp3.tsv


julia /opt/src/gatk4_m2_maflite.jl $TID $NID tmp3.tsv $PID $BUILD


filts=$(echo $DETIN_FILT | tr "," "\n")

cp $PID.m2.all.maflite.tsv  tmpA.tsv

for filt in $filts
do
	echo " [$filt]"
	grep -v $filt tmpA.tsv > tmpB.tsv
	cp tmpB.tsv  tmpA.tsv
	wc -l tmpA.tsv
done

cp tmpA.tsv $PID.m2.deTiN.maflite.tsv 




