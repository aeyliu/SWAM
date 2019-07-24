#!/usr/bin/perl -w

use strict;

#my $weights = "/net/snowwhite/home/aeyliu/pima/prediXcan/meta-TWAS-pipeline/lcl.cv.weights.dump";
my $weights = $ARGV[0];
my $genotypes_path = $ARGV[1];
my $chr = $ARGV[2];
my $out = $ARGV[3];


if(-e "$genotypes_path/chr$chr.dosage.txt.gz")
{
	my %rs2gene = ();
	my %gene2rs  = ();
	open(IN,$weights) || die "Cannot open $weights file\n";
	while(<IN>) {
	    my ($rs,$gene) = split(/\|/);
	    $rs2gene{$rs} = [] unless ( defined($rs2gene{$rs}) );
	    push(@{$rs2gene{$rs}}, $gene);
	    $gene2rs{$gene} = [] unless ( defined($gene2rs{$gene}) );
	    push(@{$gene2rs{$gene}}, $rs);    
	}
	close IN;
	
	my %rs2dose = ();
	my %genedict = ();
	my $n = 0;
	my $nind = 0;
	open(IN,"zcat $genotypes_path/chr$chr.dosage.txt.gz |") || die "Cannot open file\n";
	while(<IN>) {
	    my ($chrom,$rs,$pos,$ref,$alt,$maf,@d) = split;
	    if ( defined($rs2gene{$rs}) ) {
		$rs2dose{$rs} = \@d;
		if ( $nind == 0 ) {
		    $nind = $#d+1;
		}
		else {
		    die "$chrom $rs $pos $nind @d\n" unless ( $nind == $#d + 1 );
		}
		foreach my $g (@{$rs2gene{$rs}}) {
		    $genedict{$g} = 1;
		}
		++$n;
	    }
	}
	close IN;
	
	print STDERR "$n SNPs filtered. nind = $nind\n";
	
	open(OUT,">>$out") || die "Cannot open file\n";
	
	print OUT join("\t","GENE","RSID1","RSID2","VALUE\n");
	
	foreach my $gene (sort keys %genedict) {
	    my @rsids = ();
	    foreach my $rsid (sort @{$gene2rs{$gene}}) {
		push(@rsids,$rsid) if ( defined($rs2dose{$rsid}) );
	    }
	    print STDERR "Processing $gene across ".($#rsids+1)." rsIDs\n";    
	    for(my $i=0; $i < @rsids; ++$i) {
		my ($isum,$isq) = (0,0);
		my $ri = $rs2dose{$rsids[$i]};
		for(my $k=0; $k < $nind; ++$k) {
		    $isum += $ri->[$k];
		    #$isq += ( $ri->[$k] * $ri->[$k] );
		}
		for(my $j=0; $j <= $i; ++$j) {
		    my ($jsum,$jsq,$ijsum) = (0,0,0);		    
		    my $rj = $rs2dose{$rsids[$j]};
		    for(my $k=0; $k < $nind; ++$k) {
			$jsum += $rj->[$k];
			#$jsq += ( $rj->[$k] * $rj->[$k] );
			$ijsum += ( $ri->[$k] * $rj->[$k] );		
		    }
		    print OUT join("\t",$gene,$rsids[$j],$rsids[$i],$ijsum/$nind - $isum*$jsum/$nind/$nind)."\n";
		}
	    }
	}

	close OUT;
}