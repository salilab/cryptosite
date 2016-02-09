#!/bin/bash

awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){\
    {for(a=31;a<=38;a++){printf $a}}{printf " "}\
    {for(a=39;a<=46;a++){printf $a}}{printf " "}\
    {for(a=47;a<=54;a++){printf $a}}{printf " \n"}}' $1 |\
     awk 'BEGIN{a=0;b=0;c=0}{a+=$1;b+=$2;c+=$3}END{printf "%8.3f %8.3f %8.3f\n",a/NR,b/NR,c/NR}'
