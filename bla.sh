TT='TNi261'
ARRAY=(/scrapp/$TT/pred_dECALCrAS1000/*)
DRIN=${ARRAY[1 - 1]}
DIN=`echo ${DRIN} | cut -d '/' -f 5`
DROUT=`echo ${DRIN} | cut -d '/' -f 4,5`
DRT="./$DROUT"
PDB=`echo $DIN | cut -d \. -f 1`



echo $ARRAY
echo $DRIN
echo $DROUT
echo $DRT
echo $PDB


