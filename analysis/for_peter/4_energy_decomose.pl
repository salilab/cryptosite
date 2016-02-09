#!/usr/local/ActivePerl-5.6/bin/perl -w
#use strict

   open(INPUT,"input")|| die;
    while(<INPUT>){
             chomp;   
             @word1=split;
             $protein=$word1[0];
             $chain=$word1[1];
             $target="$protein-$chain";
             main();		   
   }
   close INPUT;


sub top{
    open (LEAP_IN, ">leapin");
               print LEAP_IN ("source /netapp/sali/AMBER/amber11/dat/leap/cmd/oldff/leaprc.ff03
#loadamberparams /netapp/sali/AMBER/amber11/dat/contrib/heme/frcmod.hemall
#loadamberprep /netapp/sali/AMBER/amber11/dat/contrib/heme/heme_all.in
set default PBradii bondi
trx = loadpdb $protein-$chain-minimize.pdb
saveamberparm trx  $target.prmtop $target\_rec.crd.1
quit\n");
   close (LEAP_IN);
   print `sleap -f leapin\n`;
}


sub sequence{ 
  open(SE,"$target-sequence")|| die"$target-sequence";
  while(<SE>){
    chomp; 
    @length=split//;
  }
  close SE;
  $length=@length;  
}


sub mmpbsa{
  open (MM, ">mm_pbsa.in");
  
  open (TEM, "mm_pbsa.in_tem");  
  while(<TEM>){
   chomp;
    if(/tem/){
      s/tem/$length/; 
      print MM "$_\n";
    }elsif(/target/){
      s/target/$target/; 
      print MM "$_\n";
    }else{
      print MM "$_\n";
    }
  }
  close TEM;
  close MM;  
}


sub main{
  
		   sequence();		   		   
		   top();
		   mmpbsa();
                   
                   print `mm_pbsa.pl  mm_pbsa.in\n`;		 		  		  	  
		   print `sander -O -i sander2.in -p $target.prmtop -c $target\_rec.crd.1  -ref $target\_rec.crd.1  -o $target.out  -r $target.res \n`; 
		   print `mv mdinfo $target\_decomp.out\n`;
		   print "------------------------------------------------------------------------------------------\n";	
}
