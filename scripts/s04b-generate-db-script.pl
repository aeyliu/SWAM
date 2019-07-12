#!/usr/bin/perl -w
use strict;

my $in = $ARGV[0];
my $out = $ARGV[1];

my $weights = "$in.weights.txt";
my $extra = "$in.extra.txt";

my $template = "/net/snowwhite/home/aeyliu/pima/prediXcan/SWAM-pipeline/utilities/sql-template.txt";
#my $extra = "/net/snowwhite/home/aeyliu/pima/prediXcan/meta-TWAS-pipeline/lcl.emp.extra.txt";
#my $weights = "/net/snowwhite/home/aeyliu/pima/prediXcan/meta-TWAS-pipeline/predicted_expression.fdr.weights.neg.txt";

open(EXTRA,"cat $extra|") || die "Cannot open file\n";

my %ensg_reg = ();
open(OUT,">>$out.sql") || die "Cannot open file\n";
open(TEMPLATE,"cat $template|") || die "Cannot open file\n";
while(<TEMPLATE>)
{
 last if(/INDEX/);
 print OUT $_;

}
close TEMPLATE;
while(<EXTRA>)
{
 #next if(/NA/);
 chomp;
 #my ($ensg,$gene,$r2,$n_snps,$pv,$qv) = split;
 my ($ensg,$gene,$n_snps) = split; #the number of entries here have to match the template file
 $ensg_reg{$ensg} = "INSERT INTO \"extra\" VALUES(\'$ensg\',\'$gene\',$n_snps)\;";
 print OUT "$ensg_reg{$ensg}\n";
}

close EXTRA;

open(WEIGHTS,"cat $weights|") || die "Cannot open file\n";
while(<WEIGHTS>)
{
 #next if(/NA/);
 chomp;
 my ($rs,$ensg,$weight,$a1,$a2) = split;
 my $line = "INSERT INTO \"weights\" VALUES(\'$rs\',\'$ensg\',$weight,\'$a1\',\'$a2\')\;";
 if(defined($ensg_reg{$ensg}))
 {
  print OUT "$line\n";
 }

}

close WEIGHTS;

open(TEMPLATE1,"cat $template|") || die "Cannot open file\n";
while(<TEMPLATE1>)
{
 next unless(/INDEX/);
 print OUT $_;

}
close TEMPLATE1;
print OUT "COMMIT\;";
close OUT;