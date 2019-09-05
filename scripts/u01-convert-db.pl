#!/usr/bin/perl -w

use strict;

my $file = $ARGV[0];
my $dir = $ARGV[1];
my $index = $ARGV[2];

#my $file = "/net/snowwhite/home/aeyliu/pima/prediXcan/git/SWAM/scripts/db-format.txt";
#my $dir = "/net/snowwhite/home/aeyliu/pima/prediXcan/git/test/v7";
#my $index = "/net/snowwhite/home/aeyliu/pima/prediXcan/git/test/v7/index.txt";

my $extra_ensg = "";
my $extra_gene = "";
my $extra_r2 = "";
my $extra_snps = "";
my $extra_pv = "";
my $extra_qv = "";
my $weights_rs = "";
my $weights_ensg = "";
my $weights_wt = "";
my $weights_ref = "";
my $weights_alt = "";

open(FORMAT,"cat $file |") || die "Cannot open file\n";
while(<FORMAT>)
{
 chomp;
 my ($type, $var, $col) = split;
 if($type eq "EXTRA")
 {
  if($var eq "ENSG") {$extra_ensg = $col}
  if($var eq "GENE") {$extra_gene = $col}
  if($var eq "R2") {$extra_r2 = $col}
  if($var eq "SNPS") {$extra_snps = $col}
  if($var eq "PV") {$extra_pv = $col}
  if($var eq "QV") {$extra_qv = $col}
 }
 if($type eq "WEIGHTS")
 {
  if($var eq "RS") {$weights_rs = $col}
  if($var eq "ENSG") {$weights_ensg = $col}
  if($var eq "WEIGHT") {$weights_wt = $col}
  if($var eq "REF") {$weights_ref = $col}
  if($var eq "ALT") {$weights_alt = $col}
 }
}

close FORMAT;

open(INDEX,"cat $index |") || die "Cannot open file\n";
while(<INDEX>)
{
 chomp;
 my @index_line = split;
 my $tissue = $index_line[0];
 open(EXTRA,"cat $dir/$tissue.extra.dump |") || die "Cannot open file\n";
 open(OUT,">> $dir/$tissue.extra.dump2") || die "Cannot open file\n";
 while(<EXTRA>)
 {
  chomp;
  my @line = split(/\|/,$_);
  print OUT join("|",$line[$extra_ensg-1],$line[$extra_gene-1],$line[$extra_r2-1],$line[$extra_snps-1],$line[$extra_pv-1],"$line[$extra_qv-1]\n");
 }
 close EXTRA;
 close OUT;

 `rm -f $dir/$tissue.extra.dump`;
 `mv $dir/$tissue.extra.dump2 $dir/$tissue.extra.dump`;

 open(WEIGHTS,"cat $dir/$tissue.weights.dump |") || die "Cannot open file\n";
 open(OUT,">> $dir/$tissue.weights.dump2") || die "Cannot open file\n";
 while(<WEIGHTS>)
 {
  chomp;
  my @line = split(/\|/,$_);
  print OUT join("|",$line[$weights_rs-1],$line[$weights_ensg-1],$line[$weights_wt-1],$line[$weights_ref-1],"$line[$weights_alt-1]\n");
 }
 close WEIGHTS;
 close OUT;

 `rm -f $dir/$tissue.weights.dump`;
 `mv $dir/$tissue.weights.dump2 $dir/$tissue.weights.dump`;
}
close INDEX;


