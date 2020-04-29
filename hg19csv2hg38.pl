#!/usr/bin/perl -w
use strict;
my $header ="clarity_sample_id,id,type,assay,sex,diagnosis,phenotype,group,father,mother,clarity_pool_id,platform,read1,read2,analysis_dir";

my $run_folder = $ARGV[0];
opendir my $dh, $run_folder or die $!;

my %csvs;
my %sample_base;
while (my $csv = readdir $dh) {
    if ($csv =~ /^\.+/) {
        next;
    }
    my $line = 0;
    my $mod_line = $header."\n";
    my $is_hg38 = 0;
    open (CSV, $run_folder.$csv) or die;
    while (<CSV>) {
        if ($line == 0) {
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

        $mod_line = $mod_line.join(',',@val[0..2]).",wgs-hg19to38,".join(',',@val[4..$#val]);

    }
    unless ($is_hg38 == 1) {
        print "new CSV  $csv\n";
        print $mod_line;
    }
}


