#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use hyunlib qw(forkExecWait);

my $directory = $ARGV[0];
my $out = $ARGV[1];

my $files = `ls $directory/*.db`;
my @files = split(/\n/, $files);


open(OUT,">$out/index.txt") || die "Cannot open $out/index.txt\n";
foreach my $file( @files )
{
 my @temp = split(/\//, $file);
 my @tissueTemp = split(/\./, $temp[$#temp]);
 pop @tissueTemp;
 my $tissue = join(".", @tissueTemp);
 print OUT join("\t", $tissue,"$file\n");
}

close OUT;