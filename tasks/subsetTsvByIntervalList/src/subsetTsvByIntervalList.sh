#!/bin/sh

BLK=""
while getopts "s:i:l:p:f:c:o:b?" Option
do
    case $Option in
        s    ) ID=$OPTARG;;
        i    ) IN=$OPTARG;;
        l    ) IL=$OPTARG;;
        c    ) CHR=$OPTARG;;
        p    ) POS=$OPTARG;;
        f    ) FST=$OPTARG;;
        b    ) BLK="-b";;
        o    ) OUT=$OPTARG;;
        ?    ) echo "Invalid option: -$OPTARG" >&2
               exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               exit 0;;
    esac
done


echo "ID:			${ID}"
echo "Intervals:	${IL}" 
echo "Input:    	${IN}" 
echo "chromsome:	${CHR}" 
echo "position:		${POS}" 
echo "stub:			${FST}" 
echo "blacklist:	${BLK}" 
echo "output area:	${OUT}" 

copt="-c $CHR "
if [[ -z $CHR ]]; then
   echo "No Chromosome field specified"
   copt=""
fi  

popt="-p $POS "
if [[ -z $POS ]]; then
   echo "No position field specified"
   popt=""
fi  

fopt="-f $FST "
if [[ -z $FST ]]; then
   echo "No stub specified"
   fopt=""
fi  


oopt="-o $OUT "
if [[ -z $OUT ]]; then
   echo "No output specified"
   oopt=""
fi   



Dir=`dirname $0`

echo ""
echo "subsetTsvByIntervalList  python command line: "
echo "python $Dir/subsetTsvByIntervalList.py -s $ID -i $IN -l $IL $copt $popt $fopt $BLK  $oopt"

python $Dir/subsetTsvByIntervalList.py -s $ID -i $IN -l $IL $copt $popt $fopt $BLK  $oopt

echo "done"
