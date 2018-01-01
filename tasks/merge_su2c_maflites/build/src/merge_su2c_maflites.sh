#!/bin/sh

BLK=""
while getopts "s:i:l:p:f:c:o:b?" Option
do
    case $Option in
        i    ) ID=$OPTARG;;
        P    ) PROD=$OPTARG;;
        2    ) M2=$OPTARG;;
        S    ) SVABA=$OPTARG;;
        s    ) STRELKA=$OPTARG;;
        t    ) TUMOR=$OPTARG;;
        n    ) NORMAL=$OPTARG;;
        ?    ) echo "Invalid option: -$OPTARG" >&2
               exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               exit 0;;
    esac
done


echo "ID:			    ${ID}"
echo "Tumor:	   ${TUMOR}" 
echo "Normal:    	  ${NORMAL}" 
echo "Production maf:	${PROD}" 
echo "M2 maflite:		${M2}" 
echo "SvABA maflite:      ${SVABA}" 
echo "Strelka maflite:      ${STRELKA}" 

Dir=`dirname $0`

echo ""
echo "merge_su2c_maflites  python command line: "
echo "python $Dir/merge_su2c_maflites.py -i $ID -t $TUMOR -n $NORMAL -P $PROD -2 $M2 -S $SVABA -s $STRELKA"

python $Dir/merge_su2c_maflites.py -i $ID -t $TUMOR -n $NORMAL -P $PROD -2 $M2 -S $SVABA -s $STRELKA

echo "done"
