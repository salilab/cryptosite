#!/usr/local/ActivePerl-5.6/bin/perl -w
#use strict

#if adding new residue for sequence, also change get_delEres.sh
@residue{"ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR","HSD","HSP","HID","HIE","HIP","CYX","HEM"}=("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","H","H","H","H","H","C","");


  $x=0;
  open(DATA,"input");
  while(<DATA>){
      chomp;
      @a=split;
      $protein=$a[0];
      $chain="A"; #chain ignored here
      $x++;
      print "$x $protein\n";
      #print `/netapp/sali/pweinkam/amber/hbplus  -h 2.7 -d 3.5 -c     $protein\n`;
      #print `/netapp/sali/pweinkam/amber/hbplus  -h 2.9 -d 4   -N -c  $protein\n`;
      print `/netapp/sali/peterc/Undrugabble/SiteCrypt/hbplus  -h 2.7 -d 3.5 -c     $protein\n`;
      print `/netapp/sali/peterc/Undrugabble/SiteCrypt/hbplus  -h 2.9 -d 4   -N -c  $protein\n`;
      sequence();

  }
  close DATA;

sub  sequence{
              open(SEQUENCE,">$protein-sequence");
             open(A,"$protein");
             while(<A>){
               chomp;
               @b=split//;
               if( defined $residue{"$b[17]$b[18]$b[19]"}){
               }else{
                    print "$_\n";
                    print "$protein not standard aminoacid\n";
                    die;
               }
               if( (/^ATOM.* CA /)){
                    $q=lc($residue{"$b[17]$b[18]$b[19]"});
                    print SEQUENCE "$q";
               }
             }
             close A;
             close SEQUENCE;

}
