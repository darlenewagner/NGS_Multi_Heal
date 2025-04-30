#!/usr/bin/perl

use strict;
use warnings;

my $idx = 0;

while(<STDIN>)
{
    $idx++;
    
    if($_ =~ /^ {8}/)
     {
	print $idx, "\n";
     }
    else
     {
	 print "\n";
     }
}

