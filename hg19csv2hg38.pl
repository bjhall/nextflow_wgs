#!/usr/bin/perl -w
use strict;
my $header ="clarity_sample_id,id,type,assay,sex,diagnosis,phenotype,group,father,mother,clarity_pool_id,platform,read1,read2,analysis_dir";

my $run_folder = $ARGV[0];
my $command = `mkdir $run_folder/pipeline_reruns19to38`;
opendir my $dh, $run_folder or die $!;

my %csvs;
my %sample_base;
while (my $csv = readdir $dh) {
    if ($csv =~ /^\.+|_/) {
        next;
    }
    my $line = 1;
    my $mod_line = $header."\n";
    my $is_hg38 = 0;
    my $sample;
    open (CSV, $run_folder.$csv) or die;
    while (<CSV>) {
        if ($line == 1) {
            $line++;
            next;
        }
        my @val = split/,/;
        if ($val[3] =~ /38/ ) {
            $is_hg38 = 1;
            next;
        }
        unless ($val[3] eq 'assay') {
            if ($sample_base{$val[0]}) {
                #print STDERR "$val[1] seen before WARN!\n";
            }
            else {
                $sample_base{$val[0]} = 1;
            }
        }
        if ($val[2] eq 'proband' && $sample_base{$val[0]}) {
            $sample = $val[1];        
        }
        my $group = '';
        my $father = '';
        my $mother = '';
        unless ($val[7] eq '') { $group = $val[7]."_38"; }
        unless ($val[8] eq '') { $father = $val[8]."_38"; }
        unless ($val[9] eq '') { $mother = $val[9]."_38"; }
        $mod_line = $mod_line.join(',',$val[0],$val[1]."_38",$val[2]).",wgs-hg19to38,".join(',',@val[4..6]).",$group,$father,$mother,".join(',',@val[10..$#val]);
        $line++;
    }
    
    unless ($is_hg38 == 1 || $sample eq '') {
        
        open (OUT, '>', $run_folder."/pipeline_reruns19to38/".$sample."_38".".csv");
        print "new CSV  $csv\n";
        print OUT $mod_line;
        close OUT;
    }
}


