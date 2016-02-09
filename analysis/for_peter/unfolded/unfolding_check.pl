#!/usr/local/ActivePerl-5.6/bin/perl -w
#use strict



#ion-->scwrl-->amber

@residue{"a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"}=("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR");




@residue=sort(keys%residue);
for($i=0;$i<=@residue-1;$i++){
       $residue=$residue[$i];   
       print "$residue{$residue}\n";
       input();	           	     		                        
}


#####################################################################################



sub input{  		  	 	
  open(LL,"$residue{$residue}/$residue{$residue}.out");
  open(CHECK,">check");
  $key1=0;$key2=0;
  while(<LL>){
    chomp;
    if(/minimize structure/){
      $key1++;
    }
          if($key1>0){
	    $key2++;
	  }
		if($key2>0){
		  print  CHECK "$_\n";
		  if($key2 == 21){
		    goto PP;
		  }
		}
  }
  PP:close LL;
  close CHECK;
  print `rm test\n`;        
  print `diff -b check sander.in_backy>test\n`;      
          $line=0;
	  $line=`wc -l < test`;
          chomp($line);	  
	  if($line != 0){
	    print "input error\n";
	    print `diff -b check sander.in_backy\n`;   
	    die;
	  }
}

