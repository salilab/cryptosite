#!/usr/local/ActivePerl-5.6/bin/perl -w
#use strict



#ion-->scwrl-->amber

@residue{"a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"}=("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR");

@residue=(keys%residue);

$AMBERHOME = "/home/.2/amber9";
		 for($i=0;$i<=@residue-1;$i++){                           
		      leapin();
		      energy_amber();  		                        		           
	          }

#####################################################################################



sub leapin{
    open (LEAP_IN, ">leapin");
               print LEAP_IN ("source leaprc.ff03
set default PBradii bondi
peptide= sequence { ACE $residue{$residue[$i]} NME}
saveamberparm peptide $residue{$residue[$i]}.top $residue{$residue[$i]}.crd
quit\n");

}


sub  energy_amber{
		      print `$AMBERHOME/exe_32/tleap -f leapin\n`;
                      print `$AMBERHOME/exe_32/sander -O -i sander.in_backy -p  $residue{$residue[$i]}.top -c $residue{$residue[$i]}.crd  -ref $residue{$residue[$i]}.crd  -o $residue{$residue[$i]}.out  -r $residue{$residue[$i]}.res \n`;
                      print `$AMBERHOME/exe_32/ambpdb -p $residue{$residue[$i]}.top <$residue{$residue[$i]}.res> $residue{$residue[$i]}.pdb\n`;                  
                      print `mv $residue{$residue[$i]}* $residue{$residue[$i]}\n`;
}

