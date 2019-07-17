#!/usr/bin/perl -w
use strict;

#this program generates the weights table for the .db aggregate predicted expression

my $correlations = $ARGV[0];
my $index_path = $ARGV[1];
my $db_path = $ARGV[2];
my $ensg_map = $ARGV[3];
my $out = $ARGV[4];

#my $index_path = "/net/snowwhite/home/aeyliu/pima/prediXcan/GTEx-V6p-1KG-2016-11-16/index.txt";
#my $correlations = "/net/snowwhite/home/aeyliu/pima/prediXcan/gtex-predictions/correlations.fdr.weighted.txt";
#my $ensg_map = "/net/snowwhite/home/aeyliu/pima/prediXcan/SWAM-pipeline/utilities/ensg.map.txt";
#my $db_path = "/net/snowwhite/home/aeyliu/pima/prediXcan/GTEx-V6p-1KG-2016-11-16";

my %emap = ();
open(ENSGMAP, "cat $ensg_map|") || die "Cannot open file\n";
while(<ENSGMAP>)
{
 chomp;
 my ($ensg,$gene) = split;
 $emap{$ensg} = $gene;
}
close ENSGMAP;


my %gene_cor_tissue = (); #hash table with correlation for each gene and each tissue (indexed by $i)
open(COR,"cat $correlations|") || die "Cannot open file\n";
while(<COR>)
{
 my ($gene,@cors) = split;
 my $i = 1;
 foreach my $cor (@cors)
 {
  $gene_cor_tissue{$gene}{$i} = $cor;
  ++$i;
 }
}


my %out_hash_table_reg = ();
my %out_hash_table_weight = ();
my %out_hash_table_a1 = ();
my %out_hash_table_a2 = ();
open(OUT,">>$out.tmp") || die "Cannot open file\n";
open(INDEX,"cat $index_path |") || die "Cannot open file\n";
my $i = 1;

while(<INDEX>)
{
 chomp;
 my @tissueLine = split;
 my $tissue = $tissueLine[0];
 print "Processing tissue $tissue\n";
 my $weights_file = "$db_path/$tissue.weights.dump";
 open(WEIGHTS,"cat $weights_file |") || die "Cannot open file\n";
 while(<WEIGHTS>)
 {
  my $line = $_;
  chomp($line);
  my ($rs,$ensg,$weight,$a1,$a2) = split(/\|/,$line);
  my $gene_temp = $emap{$ensg};
  if(defined($gene_cor_tissue{$gene_temp}{$i}))
  {
   my $cor_temp = $gene_cor_tissue{$gene_temp}{$i};
   if(defined($out_hash_table_reg{$rs}{$ensg}))
   {
    if($cor_temp>0)
    {
     $out_hash_table_weight{$rs}{$ensg} += $weight * $gene_cor_tissue{$gene_temp}{$i}; #add in the weight for the SNP multiplied by the tissue weight 
    }
    else
    {
     $out_hash_table_weight{$rs}{$ensg} += $weight * 0;
    }
   }
   unless(defined($out_hash_table_reg{$rs}{$ensg}))
   {
    $out_hash_table_reg{$rs}{$ensg} = 1;
    if($cor_temp>0)
    {
     $out_hash_table_weight{$rs}{$ensg} = $weight * $gene_cor_tissue{$gene_temp}{$i};
    }
    else
    {
     $out_hash_table_weight{$rs}{$ensg} = $weight * 0;

    }
    $out_hash_table_a1{$rs}{$ensg} = $a1;
    $out_hash_table_a2{$rs}{$ensg} = $a2;
   }
   print OUT join("\t",$rs, $ensg);
   print OUT "\n";
  }
 }
 close WEIGHTS;
 ++$i;
}
close OUT;
close INDEX;

`cat $out.tmp | sort | uniq > $out.2.tmp`;
`rm -f $out.tmp`;


my %n_snps = (); #number of SNPs per gene (for the extra table)
open(OUT,">>$out.weights.txt") || die "Cannot open file\n";
open(WTS,"cat $out.2.tmp |") || die "Cannot open file\n";
while(<WTS>)
{
 my ($rs,$ensg) = split;
 if(defined($out_hash_table_weight{$rs}{$ensg}))
 {
  unless($out_hash_table_weight{$rs}{$ensg} == 0)
  {
   $n_snps{$ensg}++;
   print OUT join("\t", $rs,$ensg,$out_hash_table_weight{$rs}{$ensg},$out_hash_table_a1{$rs}{$ensg},$out_hash_table_a2{$rs}{$ensg});
   print OUT "\n";
  }
 }
}
close OUT;
close WTS;

#make extra table 
open(OUT2,">>$out.extra.txt") || die "Cannot open file\n";
foreach my $key (keys %n_snps)
{
 print OUT2 join("\t", $key , $emap{$key}, "$n_snps{$key}\n");
}

close OUT2;


`rm -f $out.2.tmp`;
