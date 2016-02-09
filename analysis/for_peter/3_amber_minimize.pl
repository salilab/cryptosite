#!/usr/local/ActivePerl-5.6/bin/perl -w
#use strict

#ion-->scwrl-->amber
    
#$AMBERHOME = "/soft/linux/pkg/amber9-011007";
open(INPUT,"input");
   while(<INPUT>){
             chomp;            
             @word1=split;
             $protein=$word1[0];
             $chain=$word1[1];
             $o_chain=$chain;

       print "-----------------------------------------------------\n";
       print "$protein-$chain \n";      
            
          if($chain eq "XX"){
            $chain="-";
          }

          ion();#check the ionization state of HIS 
          $chain=$o_chain;
	  
            open(ORISEQ,"$protein-$chain-sequence");
                  while(<ORISEQ>){
                   chomp;
                   @sequence=split//;
                  }
            close ORISEQ;  
            $length=@sequence;    
  	      		     
		      input();		
		      leapin();
		      energy_amber();  
                               
}
close INPUT;

#####################################################################################



sub input{
                      cys();      
		      open(EF,">$protein-$chain-H_S");
		      open(EE,"$protein-$chain.pdb");			   
		      while(<EE>){
		        chomp;			
			@word77=split//;
			$tar="$word77[22]$word77[23]$word77[24]$word77[25]$word77[26]";
                        $tar=~s/ //g;			
                                 if(defined $s{$tar}){
                                     s/CYS/CYX/;
                                     print "$_\n";
                                     print EF "$_\n";
                                 }elsif(defined $ion{$tar}){
			             s/HIS/$ion{$tar}/;
			             print "$_\n";
				     print EF "$_\n";
			         }else{
			             print EF "$_\n";
			         }
		      }								
		      close EE;
		      close EF;		      		       
}




sub leapin{
    open (LEAP_IN, ">leapin");
               print LEAP_IN ("source /netapp/sali/AMBER/amber11/dat/leap/cmd/oldff/leaprc.ff03
#loadamberparams /netapp/sali/AMBER/amber11/dat/contrib/heme/frcmod.hemall
#loadamberprep /netapp/sali/AMBER/amber11/dat/contrib/heme/heme_all.in
set default PBradii bondi
trx = loadpdb $protein-$chain-H_S
saveamberparm trx $protein-$chain.top $protein-$chain.crd
quit\n");	       	       	     
   close (LEAP_IN);
}




sub  energy_amber{		      
		      print `sleap -f leapin\n`;	  
                      print `sander -O -i sander.in -p $protein-$chain.top -c $protein-$chain.crd  -ref $protein-$chain.crd  -o $protein-$chain.out  -r $protein-$chain.res \n`;
                      print `ambpdb -p $protein-$chain.top <$protein-$chain.res> $protein-$chain-minimize.pdb\n`;
		      print `cp $protein-$chain.out $protein-$chain-minimize.out\n`;

}

sub ion{#analysis the output filr generated from HBPLUS
  open(DD,"$protein-$chain.hb2"); 
  $aa="  A "; $tt="  T "; $cc="  C "; $gg="  G "; 
  while(<DD>){
    chomp;#A0017-HIS ND1 Z0002-HOH O
    if(/$aa/){     
      goto BY;
    }elsif(/$tt/){    
      goto BY;
    }elsif(/$cc/){
      goto BY;
    }elsif(/$gg/){
      goto BY;
    }
    @word22=split;
    if($word22[0]=~/.[0-9][0-9][0-9][0-9]/){#caring about the donor parts of HIS 
      if($word22[0]=~/HIS/){
         if("$word22[1]" ne "N"){#caring about the ND1 and NE2 parts of HIS 
	    print "$protein $_ \n";	    	   	       
	   @word55=split;#$word55[0]="A0017-HIS"
           if( ($word55[2]=~/...../) || ($word55[2]=~/HOH/)){
	      print "$protein $_ \n";
    
	      @word44=split(//,$word55[0]);#A 0 0 0 0 - H I S
	      ($a,@word44)=@word44;	   #  0 0 0 0 - H I S


	   #   while($word44[0] =~/[^1-9]/){
	   #     ($a,@word44)=@word44;	       
	   #   }#@word44=["1","7","-","H","I,"S"];  		       

                     for($gh=1;$gh<=3;$gh++){
                         if($word44[0] =~/[^1-9]/){
                          ($a,@word44)=@word44;
                         }
                     }

		   print "@word44\n";		   
		   for($pp=0;$pp<=2;$pp++){
		    $a=pop(@word44);
		   }
                   if($word44[@word44-1] eq "\-"){
                       $a=pop(@word44);
                   }
		   print "@word44\n";		   
		   $target=join("",@word44);#17
		   $ion{$target}.=$word55[1];#$word55[1]="ND1"
	           print "$target $ion{$target}\n";
	 	 
	 }	 	 
	 }
      }
    } 
      BY:;
  }
  close DD; 
		     #$ion{137}="ND1ND1";
		     print "********************************************************************\n";
		     @key=sort{$a<=>$b}keys(%ion);
		     
		     foreach(@key){
		       $first=$_;
		       $last=$ion{$first};
		     
		       if($last=~/ND1/){
		          if($last=~/NE2/){
			    $ion="HIP";
			    $ion{$first}=$ion;
			    
			    print "$protein $first $ion{$first}\n";
			  }else{
			    $ion="HID";
			    $ion{$first}=$ion;
			    print "$protein $first $ion{$first}\n";
			  }

		       }else{
		            $ion="HIE"; 
			    $ion{$first}=$ion;
		            print "$protein $first $ion{$first}\n";
		       }
		     }
		      print "*******************************************************************\n";
 
}


sub cys{
  undef %s;
  open(S,"Sulfate");
  while(<S>){
     chomp;
     @words=split;
     if(/$protein-./){
         $s{$words[0]}++;
     }
  }
  close S;
}
