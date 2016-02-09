#!/usr/bin/perl -w

use strict;
use List::Util 'shuffle';

my @lines = <>;
print shuffle( @lines ); 
