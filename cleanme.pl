#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $FULL_CLEAN = $ARGV[0] || 0;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";




my @files_to_keep = qw (
README.md
runTrinityDemo.pl
cleanme.pl
);                      


my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (! $keep{$file}) {
        if (-f $file) {
            print STDERR "-removing file: $file\n";
            unlink($file);
        }
        
    }
}

`rm -rf ./wt_*`;
`rm -rf ./GSNO_*`;
`rm -rf ./trinity_out_dir`;
`rm -rf ./edgeR`;
`rm -rf ./__TrinDemo_checkpoints_dir`;

`rm -f ./data/*fastq`;
`rm -f ./data/*readcount`;

exit(0);
