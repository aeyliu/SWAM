#!/usr/bin/perl -w

use strict;

my $indexf = $ARGV[0];
my $sql_path = $ARGV[1];
my $out_path = $ARGV[2];

open(INDEX,"cat $indexf |") || die "Cannot open index file!\n";
while(<INDEX>)
{
 chomp;
 my ($name,$file) = split;
 `$sql_path $file 'SELECT rsid, gene, weight, ref_allele, eff_allele FROM weights' > $out_path/$name.weights.dump`;
 `$sql_path $file 'SELECT gene, genename, [pred.perf.R2], [n.snps.in.model], [pred.perf.pval] FROM extra' > $out_path/$name.extra.dump`;

}

close INDEX;