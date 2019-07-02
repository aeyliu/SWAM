#!/usr/bin/perl -w

use strict;
use Cwd qw(cwd);

my $outPrefix = "/net/snowwhite/home/aeyliu/pima";

print "$outPrefix\n";

my @outArray = split(/\//,$outPrefix);
my $outName = $outArray[$#outArray];
my $outDir = "";
print "$#outArray\n";
if($#outArray > 0)
{
 $outDir = join("/", @outArray[0 .. $#outArray-1]);
}

print "$outName\n";
print "$outDir\n";

my $wd = cwd;

if($outDir eq "")
{
 print join("/", $wd, "$outName\n");
}

unless($outDir eq "")
{
 print join("/",$outDir,"$outName\n");
}

