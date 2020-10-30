#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;
use List::MoreUtils qw(first_index);

my %opt = ();
GetOptions( \%opt, 'vcf=s', 'ped=s', 'yaml=s'  );

my $pedfile = $opt{ped};
my ($PED, $proband, $pedsize) = read_ped($pedfile);
my $father;
my $mother;
if ($pedsize > 2) {
	$mother = $PED->{$proband}->{MOTHER};
	$father = $PED->{$proband}->{FATHER};
}
#print Dumper($PED);
my $vcf = CMD::vcf2->new('file'=>$opt{vcf} );

my @header = split/\n/,$vcf->{header_str};
my @sample_order = (sort keys %{ $PED});
foreach my $header (@header) {
	next if ($header =~ /^##INFO=<ID=(Compounds|GeneticModels|ModelScore|RankScore|RankResult)/);
	next if ($header =~ /^##Software=<ID=genmod/);
	if ($header =~ /^#CHROM/) {
		print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
		print join("\t", @sample_order);
		print "\n";
	}
	else {
		print $header."\n";
	}
	
}

my $rankresult_meta = $vcf->{meta}->{INFO}->{RankResult}->{Description};
#print $rankresult_meta."\n";
$rankresult_meta =~ s/\|/:/g;


while ( my $a = $vcf->next_var() ) {
	my ($vcfedited,$ok) = vcfstr($a, \@sample_order);
	if ($ok) {
		print $vcfedited;
	}
	else {
		next;
	}

}

sub read_ped {
	my $pedfile = shift;

	open (PED, $pedfile) or die $!;
	
	my $proband;
	my $mother;
	my $father;
	my %PED;
	
	while ( <PED> ) {
		my @line = split/\t/,$_;
		my %ind;
		
		$ind{FATHER} = $line[2];
		$ind{MOTHER} = $line[3];
		$ind{SEX} = $line[4];
		$ind{PHENO} = $line[5];
		$PED{$line[1]} = \%ind;

		unless ($line[2] eq "0" && $line[3] eq "0") {
			$proband = $line[1];
			$father = $line[2];
			$mother = $line[3];
		}
	}
	my $count = keys %PED;
	## if single sample, proband is obvious.
	if ($count <= 2 ) {
		foreach my $ind (keys %PED) {
			$proband = $ind;
		}
	}
	return \%PED, $proband, $count;
}

sub vcfstr {
	my( $v, $sample_order ) = @_;
	
	my @all_info;
	my $tot_str = $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

	# Generate and print INFO field AND STRIP SCORES
	for my $info_key (@{$v->{INFO_order}}) {
		next if ($info_key eq 'Compounds');
		next if ($info_key eq 'GeneticModels');
		next if ($info_key eq 'ModelScore');
		next if ($info_key eq 'RankScore');
		next if ($info_key eq 'RankResult');
		if($info_key eq "CSQ") {
			push @all_info, $info_key."=".$v->{_CSQstr};
		}
		else {
			push @all_info, $info_key."=".$v->{INFO}->{$info_key};
		}
	}
	$tot_str = $tot_str.join(";", @all_info)."\t";

	# Print FORMAT field
	$tot_str = $tot_str.join(":", @{$v->{FORMAT}})."\t";


	my %order;
	my $i=1;
	$order{$_} = $i++ foreach @{$sample_order};
	# Only include individuals in PED
	my @GTS;
	foreach my $sample (@{$v->{GT}}) {
		if ($order{$sample->{_sample_id}}) {
			push @GTS,$sample;
		}
	}
	# Print GT fields for all samples
	my $onevariantplease = 0;
	for my $gt ( sort {$order{$a->{_sample_id}} <=> $order{$b->{_sample_id}}} @GTS) {
		my @all_gt;
		
		for my $key ( @{$v->{FORMAT}} ) {
			## scout wont load GT-fields = "." Need to set dummy-value
			if ($key eq 'GT') {
				if ($gt->{$key} =~ /(0\/1|\.\/1|1\/1)/) {
					$onevariantplease = 1;
				}
			}

			if ($gt->{$key} eq '.') {
				push @all_gt,"0,0";
			}
			else {
				push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
			}
			#push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");

		}
		$tot_str = $tot_str.join(":", @all_gt)."\t";
	}
	$tot_str = $tot_str."\n";
	return $tot_str,$onevariantplease;
}
