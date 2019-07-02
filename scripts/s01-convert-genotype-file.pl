#!/usr/bin/perl -w
use strict;
use FindBin;
use wGetOptions
use List::Util qw(sum);

my $genotypes_file = $ARGV[0];
my $output_directory = $ARGV[1];
#my $genotypes_file = "/net/snowwhite/home/aeyliu/pima/prediXcan/gtex-predictions/snps/Cells_EBV-transformed_lymphocytes_Analysis.snps.txt";
my $dbsnp = "/net/snowwhite/home/aeyliu/pima/prediXcan/SWAM-pipeline/utilities/dbsnp_142.b37.rssorted.txt.gz";
#my $output_directory = "/net/snowwhite/home/aeyliu/pima/prediXcan/gtex-predictions/genotypes";

my %rsids = ();
print "Processing DBSNP...\n";
open(DBSNP, "zcat $dbsnp|") || die "Cannot open file\n";
while(<DBSNP>)
{
 next unless (/b142/);
 my ($build, $rs, $chr, $pos, $ref, $alt, $info) = split;
 $rsids{$chr}{$pos} = $rs;
}
close DBSNP;

print "Processing genotype file...\n";
my $current_chr = 1;
open(GENO, "cat $genotypes_file|") || die "Cannot open file\n";
open(OUT, ">>$output_directory/chr$current_chr.dosage.txt") || die "Cannot open file\n";
open(OUT1,">>$output_directory/samples.tmp") || die "Cannot open file\n";

my $line = <GENO>;
my ($id,@samples) = split("\t",$line);
print OUT1 join("\n", @samples);
close OUT1;

while(<GENO>)
{
 next if(/Id/);

 my ($variant, @dosages) = split;
 my ($chr, $pos, $ref, $alt, $build) = split("_",$variant);
 if($chr != $current_chr)
 {
  close OUT;
  $current_chr = $chr;
  open(OUT, ">>$output_directory/chr$current_chr.dosage.txt") || die "Cannot open file\n";
 }
 if(defined($rsids{$chr}{$pos}))
 {
  print OUT join("\t", "chr$chr", "rs$rsids{$chr}{$pos}", $pos, $ref, $alt, sum(@dosages)/($#dosages*2),@dosages);
  print OUT "\n";
 }
}

close GENO;
close OUT;

print "Zipping files...\n";
my @chrs = (1..22);
foreach my $chr (@chrs)
{
 `bgzip $output_directory/chr$chr.dosage.txt`;
}


`paste $output_directory/samples.tmp $output_directory/samples.tmp > $output_directory/samples.txt`;

