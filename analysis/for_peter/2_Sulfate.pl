#!/usr/local/ActivePerl-5.6/bin/perl -w
#use strict

  $x=0;
  open(DATA,"input");

  while(<DATA>){
      chomp;
      @word1=split;
      $protein=$word1[0];
      $chain=$word1[1];
        
      $x++;
      print "$x\n";
      open(S,">Sulfate");
 
      print `grep -h \"CSS.*CSS\" $protein-$chain.nb2>test\n`;     
      open(TEST,"test");
      while(<TEST>){
            chomp;
            @word2=split;
            $a=$word2[0];           $b=$word2[2];
            @worda=split(//,$a);    @wordb=split(//,$b);
            $a="$worda[1]$worda[2]$worda[3]$worda[4]$worda[5]";
            $a=~s/\-//g;

            while($a=~/\b0/){
                $a=~s/\b0//g;  
            } 

            if($worda[0] eq "-"){
                $worda[0]="XX";
            }
            $target{"$a $protein-$worda[0] $a"}++; 

            $b="$wordb[1]$wordb[2]$wordb[3]$wordb[4]$wordb[5]";
            $b=~s/\-//g;

            while($b=~/\b0/){
                $b=~s/\b0//g;
            }

            if($wordb[0] eq "-"){
                $wordb[0]="XX";
            }      
            $target{"$b $protein-$wordb[0] $b"}++;
            print "$_\n";
      }
      close TEST;  
 
            @target=sort(keys%target);
            for($i=0;$i<=@target-1;$i++){
              print S "$target[$i]\n";
            }
            undef %target;
            close S;
 }

 close DATA;
