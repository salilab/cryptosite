#!/usr/local/ActivePerl-5.6/bin/perl -w
#use strict

@residue{"ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"}=("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");

  $x=0;
  open(DATA,"input");
  while(<DATA>){
      chomp;
      @a=split;
      $protein=$a[0];$chain=$a[1];
      $x++;
      print "$x $protein\n";
      get_single_chain(); 
      last_line_check();
      aminoacid_check();
  }
  close DATA;

sub get_single_chain{
     open(P,">$protein-$chain.pdb");
     open(PDB,"pdb$protein.ent") || die "pdb$protein.ent\n";
     while(<PDB>){
      chomp;
      @b=split//;
      @bb=split;

      if(/^ATOM/){
      if( ( ($b[16] eq " ") || ($b[16] eq "A") || ($b[16] eq "1")) && ($bb[2]!~/^H/) ){
        if( $b[21] eq "$chain"){
          print P "$_\n";
        }
      }
     }
    }
     close PDB;
     close P;
}


sub last_line_check{
             open(A,"$protein-$chain.pdb");
             while(<A>){
               chomp;
               $target=$_;
             }#the last line
             close A;
             if($target=~/[0-9]  N /){
                    print "PDB last line error\n";
             }
             close A;
}

sub aminoacid_check{
             open(SEQUENCE,">$protein-$chain-sequence");
             open(A,"$protein-$chain.pdb");
             while(<A>){
               chomp;
               @b=split//;
               if( defined $residue{"$b[17]$b[18]$b[19]"}){
               }else{
                    print "$_\n";
                    print "$protein-$chain not standard aminoacid\n";
                    die;
               }
               if( (/^ATOM.* CA /) && ($b[21] eq "$chain")){
                    $q=lc($residue{"$b[17]$b[18]$b[19]"});
                    print SEQUENCE "$q";
               }
             }
             close A;
             close SEQUENCE;
}
