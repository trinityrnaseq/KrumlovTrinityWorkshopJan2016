#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use FindBin;
use Cwd;

use lib ("$FindBin::Bin/PerlLib");
use IniReader;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

#################################################################################
#
#  --autopilot         automatically run the pipeline end-to-end
#
#################################################################################


__EOUSAGE__

    ;


my $help_flag = 0;
my $AUTO_MODE = 0;


&GetOptions( 'help|h' => \$help_flag,
             'autopilot' => \$AUTO_MODE,
             
    );

if ($help_flag) {
    die $usage;
}



my $trinity_dir = $ENV{TRINITY_HOME} or die "Error, need env var TRINITY_HOME set to Trinity installation directory";
$ENV{PATH} .= ":$trinity_dir";  ## adding it to our PATH setting.

my $trinotate_dir = $ENV{TRINOTATE_HOME} or die "Error, need env var TRINOTATE_HOME set to Trinotate installation directory (note this is different than Trinity) ";


my $OS_type = `uname`;

## first check for tools needed.

my @tools = qw (Trinity
    bowtie
    samtools
    igv.sh
);
 
{
    my $missing_tool_flag = 0;
    foreach my $tool (@tools) {
        my $path = `which $tool`;
        unless ($path =~ /\w/) {
            print STDERR "Error, cannot find path to: $tool\n";
            $missing_tool_flag = 1;
        }
    }
    
    if ($missing_tool_flag) {
        die "\n\nTools must be in PATH setting before proceeding.\n\n";
    }

}



{ 
    ## unzip the gzipped files.
    foreach my $file (<data/*.gz>) {
        my $unzipped_file = $file;
        $unzipped_file =~ s/\.gz//;
        unless (-s $unzipped_file) {
            my $ret = system("gunzip -c $file > $unzipped_file");
            if ($ret) {
                die "Error, could not gunzip file $file";
            }
        }
    }
}

my $checkpoints_dir = $FindBin::Bin . "/__TrinDemo_checkpoints_dir";
unless (-d $checkpoints_dir) {
    mkdir $checkpoints_dir or die "Error, cannot mkdir $checkpoints_dir";
}


my %samples = ( 
    'GSNO' => {
        'GSNO_SRR1582646' => [ 'data/GSNO_SRR1582646_1.fastq',
                               'data/GSNO_SRR1582646_2.fastq' ],
        'GSNO_SRR1582647' => [ 'data/GSNO_SRR1582647_1.fastq',
                               'data/GSNO_SRR1582647_2.fastq' ],
        'GSNO_SRR1582648' => [ 'data/GSNO_SRR1582648_1.fastq',
                               'data/GSNO_SRR1582648_2.fastq' ],
    },
    
    'WT' => {
        'wt_SRR1582649' => [ 'data/wt_SRR1582649_1.fastq',
                             'data/wt_SRR1582649_2.fastq' ],
        'wt_SRR1582650' => [ 'data/wt_SRR1582650_1.fastq',
                             'data/wt_SRR1582650_2.fastq' ],
        'wt_SRR1582651' => [ 'data/wt_SRR1582651_1.fastq',
                             'data/wt_SRR1582651_2.fastq' ]
    } 
    
    );


##############
# run Trinity.

my ($left_fq_files_aref, $right_fq_files_aref) = &get_fq_files_listings(%samples);

my $run_Trinity_cmd = "$trinity_dir/Trinity --seqType fq "
    . " --left " . join(",", @$left_fq_files_aref)
    . " --right " . join(",", @$right_fq_files_aref) 
    . " --CPU 4 --max_memory 2G --min_contig_length 150 "; 


&process_cmd($run_Trinity_cmd, "$checkpoints_dir/trinity.ok");

# count the number of transcripts assembled.
&process_cmd("grep '>' trinity_out_dir/Trinity.fasta | wc -l ", "$checkpoints_dir/count_trans.ok");

# Examine top of Trinity.fasta file
&process_cmd("head trinity_out_dir/Trinity.fasta", "$checkpoints_dir/head_trinity.ok");

# Get Trinity stats:
&process_cmd("$trinity_dir/util/TrinityStats.pl trinity_out_dir/Trinity.fasta", "$checkpoints_dir/trin_stats.ok");


###################################
## Abundance estimation using RSEM
###################################

# also, write the samples.txt file

open (my $ofh, ">samples.txt") or die "Error, cannot write to file samples.txt";
my @rsem_result_files;
my @bam_files;

foreach my $condition (sort keys %samples) {
    
    my $samples_href = $samples{$condition};
    
    foreach my $sample (keys %$samples_href) {
        
        print $ofh join("\t", $condition, $sample) . "\n";
        
        my $replicates_aref = $samples_href->{$sample};
        
        my ($left_fq, $right_fq) = @{$replicates_aref};
        
        my $output_dir = "$sample.RSEM";
        
        my $rsem_result_file = "$output_dir/$sample.isoforms.results";
        push (@rsem_result_files, $rsem_result_file);
        
        
        my $align_estimate_command = "$trinity_dir/util/align_and_estimate_abundance.pl --seqType fq "
            . " --left $left_fq --right $right_fq "
            . " --transcripts trinity_out_dir/Trinity.fasta "
            . " --output_prefix $sample --est_method RSEM "
            . " --aln_method bowtie --trinity_mode --prep_reference --coordsort_bam "
            . " --output_dir $output_dir";
        
        &process_cmd($align_estimate_command, "$checkpoints_dir/$sample.align_estimate.ok");
        
        push (@bam_files, "$output_dir/$sample.bowtie.csorted.bam");
        
        # look at the output
        &process_cmd("head $rsem_result_file", "$checkpoints_dir/head.$sample.rsem.ok");
    }
    
    
}
close $ofh; # samples.txt

## generate matrix of counts and perform TMM normalization
&process_cmd("$trinity_dir/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix Trinity_trans @rsem_result_files", "$checkpoints_dir/counts_matrix.ok");

## Look at the counts matrix
&process_cmd("head -n20 Trinity_trans.counts.matrix", "$checkpoints_dir/head.counts.matrix.ok");

## Look at the expression matrix:
&process_cmd("head -n20 Trinity_trans.TMM.EXPR.matrix", "$checkpoints_dir/head.expr.matrix.ok");

## Examine the E90N50 statistic
&process_cmd("$trinity_dir//util/misc/contig_ExN50_statistic.pl  Trinity_trans.TMM.EXPR.matrix trinity_out_dir/Trinity.fasta > ExN50.stats", "$checkpoints_dir/ExNstats.ok");

&process_cmd("cat ExN50.stats", "$checkpoints_dir/cat_ExNstats.ok");

## Plot the values
&process_cmd("$trinity_dir/util/misc/plot_ExN50_statistic.Rscript ExN50.stats", "$checkpoints_dir/plot_ExN50.ok");

&show("ExN50.stats.plot.pdf");

## Examine read alignments in IGV
&process_cmd("igv.sh -g trinity_out_dir/Trinity.fasta " . join(",", @bam_files), "$checkpoints_dir/igv_trinity_reads.ok");



##############
## DE analysis
##############

## make a samples file
&process_cmd("cat samples.txt", "$checkpoints_dir/examine_samples_txt.ok");

## run edgeR
&process_cmd("$trinity_dir/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Trinity_trans.counts.matrix --samples_file samples.txt --method edgeR --output edgeR", "$checkpoints_dir/run.edgeR.ok");

# take a look at what edgeR generated:
&process_cmd("ls -ltr edgeR/", "$checkpoints_dir/ls.edgeR.dir.ok");


&process_cmd("head edgeR/Trinity_trans.counts.matrix.GSNO_vs_WT.edgeR.DE_results", "$checkpoints_dir/head.edgeR.DE_results.ok");

&show("edgeR/Trinity_trans.counts.matrix.GSNO_vs_WT.edgeR.DE_results.MA_n_Volcano.pdf");


eval {
    &process_cmd("cd edgeR", "$checkpoints_dir/cd.edgeR.ok");
};
if ($@) {
    system("touch $checkpoints_dir/cd.edgeR.ok");
}

# now do it in the script. :)
chdir("edgeR") or die "Error, could not cd to edgeR/"; 
print STDERR "\n ** Note: if you see an error message above about not being able to cd, just ignore it... it's a weird behavior of this demo script. Rest assured we've 'cd edgeR' just fine.   :)\n\n";

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt -P 1e-3 -C 2",
             "$checkpoints_dir/analyze_diff_expr.ok");


&process_cmd("wc -l diffExpr.P1e-3_C2.matrix", "$checkpoints_dir/wc_diff_expr_matrix.ok"); # number of DE transcripts + 1

&show("diffExpr.P1e-3_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData",
             "$checkpoints_dir/cut_clusters_tree.ok");
&show("diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60/my_cluster_plots.pdf");

print STDERR "\n\n\tDemo complete.  Congratulations!  :)\n\n\n\n";


exit(0);

####
sub process_cmd {
    my ($cmd, $checkpoint) = @_;

    unless ($checkpoint) {
        die "Error, need checkpoint file defined";
    }
    
    if (-e $checkpoint) { return; }

    
    unless ($AUTO_MODE) {
        
        my $response = "";
        while ($response !~ /^[YN]/i) {
            print STDERR "\n\n"
                . "###############################################\n"
                . "CMD: $cmd\n"
                . "###############################################\n\n"
                . "Execute (Y/N)? ";

            $response = <STDIN>;
        }

        if ($response =~ /^N/i) {
            print STDERR "\t *** Exiting on demand. ****\n\n"
                . "Goodbye. \n\n\n";
            exit(0);
        }
    }
    
    print STDERR "\tNow running:\n\t\t$cmd\n\n\n";
    
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    system("touch $checkpoint");
    
    return;
}


sub show {
    my ($image) = @_;

    my $cmd;

    if ($OS_type =~ /linux/i) {
        ## use xpdf
        $cmd = "xpdf $image";
    }
    else {
        ## rely on ImageMagick:
        $cmd = "open $image";
    }
    
    if ($AUTO_MODE) {
        $cmd .= " & ";
    }
    
    &process_cmd($cmd, "$checkpoints_dir/view." . basename($image) . ".ok");

    return;
}



####
sub get_fq_files_listings {
    my (%samples) = @_;
    
    my @left_fq_files;
    my @right_fq_files;

    foreach my $condition_fq_lists_href (values %samples) {
        
        foreach my $replicates_fq_lists_aref (values %$condition_fq_lists_href) {
            
            my ($left_fq, $right_fq) = @$replicates_fq_lists_aref;
            push (@left_fq_files, "$left_fq");
            push (@right_fq_files, "$right_fq");
        }
    }

    return(\@left_fq_files, \@right_fq_files);
}


