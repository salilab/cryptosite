#!/bin/bash

FIL_MMPBSA=$1
FIL_PDB=$2 #to get sequence, assumed to contain only chain/residues 
FIL_FREE=/netapp/sali/pweinkam/amber/unfolded/EFREE

#load energies for free amino acid
INT=(null `awk '{print $7}' $FIL_FREE`)
VDW=(null `awk '{print $13}' $FIL_FREE`)
ELE=(null `awk '{print $19}' $FIL_FREE`)
GAS=(null `awk '{print $25}' $FIL_FREE`)
GB=(null `awk '{print $31}' $FIL_FREE`)
GBSUR=(null `awk '{print $37}' $FIL_FREE`)
GBSOL=(null `awk '{print $43}' $FIL_FREE`)
GBTOT=(null `awk '{print $49}' $FIL_FREE`)

#load sequence
SEQ=(`awk '($3=="CA"){print $4}' $FIL_PDB | awk '\
($1=="ALA"){print "1"}($1=="ARG"){print "2"}($1=="ASN"){print "3"}($1=="ASP"){print "4"}\
($1=="CYS"){print "5"}($1=="CYX"){print "5"}($1=="GLN"){print "6"}($1=="GLU"){print "7"}($1=="GLY"){print "8"}\
($1=="HIS"||$1=="HID"||$1=="HIE"||$1=="HIP"||$1=="HSD"||$1=="HSP"){print "9"}\
($1=="ILE"){print "10"}($1=="LEU"){print "11"}($1=="LYS"){print "12"}\
($1=="MET"){print "13"}($1=="PHE"){print "14"}($1=="PRO"){print "15"}($1=="SER"){print "16"}\
($1=="THR"){print "17"}($1=="TRP"){print "18"}($1=="TYR"){print "19"}($1=="VAL"){print "20"}'`)

sINT=(`awk '(NR>4){print $7}' $FIL_MMPBSA`)
sVDW=(`awk '(NR>4){print $13}' $FIL_MMPBSA`)
sELE=(`awk '(NR>4){print $19}' $FIL_MMPBSA`)
sGAS=(`awk '(NR>4){print $25}' $FIL_MMPBSA`)
sGB=(`awk '(NR>4){print $31}' $FIL_MMPBSA`)
sGBSUR=(`awk '(NR>4){print $37}' $FIL_MMPBSA`)
sGBSOL=(`awk '(NR>4){print $43}' $FIL_MMPBSA`)
sGBTOT=(`awk '(NR>4){print $49}' $FIL_MMPBSA`)

RANGES=`count.pl 0 $((${#SEQ[@]} - 1))`
for r in $RANGES; do
    res=${SEQ[$r]}
    #INT VDW ELE GBSOL GBTOT
    echo $(($r + 1)) `pbc "${sINT[$r]}-(${INT[$res]})"` `pbc "${sVDW[$r]}-(${VDW[$res]})"` `pbc "${sELE[$r]}-(${ELE[$res]})"` `pbc "${sGB[$r]}-(${GB[$res]})"` `pbc "${sGBSUR[$r]}-(${GBSUR[$res]})"` `pbc "${sGBTOT[$r]}-(${GBTOT[$res]})"` | awk '{printf "%4.0f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f\n",$1,$2,$3,$4,$5,$6,$7}'
done

