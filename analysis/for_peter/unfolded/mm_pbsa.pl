#!/usr/bin/perl -w
#
# Script to perform mm_pbsa
#
# Following the idea of Irina Massova,
#   however 
#     * parameter input is now done via an input file
#     * MM energies are calc with sander
#     * GB energies are calc with sander
#     * SAS may be calc with molsurf or sander
#     * gasphase entropies are calc with nmode
#     * energy decomposition into residue-related and 
#         pairwise residue-related contributions
#         is possible
#
# Holger Gohlke
#   Last change: 16.01.2002
#
########################################################################

use strict;
use lib ("$ENV{AMBERHOME}/src/mm_pbsa");
use mm_pbsa_global       qw(@GC_PAR);
use mm_pbsa_util         qw(init_data final_clean_up);
use mm_pbsa_readinput    qw(read_input check_sanity);
use mm_pbsa_createinput  qw(create_input);
use mm_pbsa_createcoords qw(create_coords);
use mm_pbsa_calceneent   qw(calc_energy_entropy);
use mm_pbsa_statistics   qw(calc_stat);

########################################################################

MAIN:{
######

  $ | = 1; # Autoflush

  if(scalar(@ARGV) != 1){
    print "\n";
    print "USAGE: mm_pbsa.pl <input file>\n";
    die("\n");
  }
 
  # Variable declarations
  my (%GEN,%DEL,%DEC,%SAN,%GBO,%MOL,%NMO,%MAK,%PRO);
  my $TRA = "";

  # Initialize some data
  &init_data(\%PRO);

  # Read input data
  &read_input($ARGV[0],\%GEN,\%DEC,\%SAN,\%DEL,\%GBO,\%MOL,\%NMO,\%MAK,\%PRO,\$TRA);

  # Check sanity
  &check_sanity(\%GEN,\%DEC,\%SAN,\%DEL,\%GBO,\%MOL,\%NMO,\%MAK,\%PRO,\$TRA);

  # Create input files for programs
  &create_input(\%GEN,\%DEC,\%SAN,\%DEL,\%GBO,\%MOL,\%NMO,\%MAK,\%PRO,\$TRA);

  if($GEN{"GC"} || $GEN{"AS"}){
    # Creating coordinates out of trajectory
    &create_coords(\%PRO,\$TRA);
  }
  
  if($GEN{"MM"} || $GEN{"NM"} || $GEN{"PB"} || $GEN{"GB"} || $GEN{"MS"}){
    # Calc single energy and entropy contributions
    &calc_energy_entropy(\%GEN,\%DEC,\%DEL,\%MOL,\%PRO);
  }

  if($GEN{"MM"} || $GEN{"PB"} || $GEN{"GB"} || $GEN{"MS"} || $GEN{"NM"}){
    # Calc statistics
    &calc_stat(\%GEN,\%DEC,\%DEL,\%GBO,\%PRO);
  }

  # Final clean up
  &final_clean_up(\%GEN,\%DEL);
}
