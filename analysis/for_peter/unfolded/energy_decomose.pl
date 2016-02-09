#!/usr/local/ActivePerl-5.6/bin/perl -w
#use strict

@residue{"a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"}=("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR");
@residue=sort(keys%residue);

   
   for($i=0;$i<=@residue-1;$i++){                                     
	     $target=$residue{$residue[$i]};
	     print "$target\n";
             main();		   
   }
   


sub top{
    open (LEAP_IN, ">leapin");
               print LEAP_IN ("source leaprc.ff03
set default PBradii bondi
trx = loadpdb $target/$target.pdb
saveamberparm trx  $target.prmtop $target\_rec.crd.1
quit\n");
   close (LEAP_IN);
   print `/home/.2/amber9/exe_32/tleap -f leapin\n`;
}


sub sequence{
  $length=3;  
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
                   print `./mm_pbsa.pl*  mm_pbsa.in\n`;		 		  		  	  
		   print `/home/.2/amber9/exe_32/sander -O -i backy -p $target.prmtop -c $target\_rec.crd.1  -ref $target\_rec.crd.1  -o $target.de.out  -r $target.res \n`;   	          		   		   
		   print "------------------------------------------------------------------------------------------\n";	
                   print `mv $target* $target\n`;
                   print `mv mdinfo $target\n`;

}
