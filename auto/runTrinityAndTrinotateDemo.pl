#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use FindBin;
use Cwd;

use lib ("$FindBin::Bin/PerlLib");
use IniReader;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $trinotate_conf_file = "$FindBin::Bin/conf.txt";

my $usage = <<__EOUSAGE__;

#################################################################################
#
#  Optional:
#
#  --autopilot         automatically run the pipeline end-to-end
#
#  --conf <string>     configuration file to use.  Default: $trinotate_conf_file
#
#################################################################################


__EOUSAGE__

    ;


my $help_flag = 0;
my $AUTO_MODE = 0;


&GetOptions( 'help|h' => \$help_flag,
             'autopilot' => \$AUTO_MODE,
             
             'conf=s' => \$trinotate_conf_file,
             
    );

if ($help_flag) {
    die $usage;
}


unless ($trinotate_conf_file =~ m|^/|) {
    $trinotate_conf_file = cwd() . "/$trinotate_conf_file";
}

my $trinity_dir = $ENV{TRINITY_HOME} or die "Error, need env var TRINITY_HOME set to Trinity installation directory";
$ENV{PATH} .= ":$trinity_dir";  ## adding it to our PATH setting.

my $trinotate_dir = $ENV{TRINOTATE_HOME} or die "Error, need env var TRINOTATE_HOME set to Trinotate installation directory (note this is different than Trinity) ";



my $OS_type = `uname`;

my $workdir = cwd();

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


my $checkpoints_dir = $workdir . "/__TrinDemo_checkpoints_dir";
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


my $STEP_COUNT = 0; # incremented at each process_cmd


##############
# run Trinity.

my ($left_fq_files_aref, $right_fq_files_aref) = &get_fq_files_listings(%samples);

my $run_Trinity_cmd = "$trinity_dir/Trinity --seqType fq "
    . " --left " . join(",", @$left_fq_files_aref)
    . " --right " . join(",", @$right_fq_files_aref) 
    . " --CPU 4 --max_memory 2G --min_contig_length 150 "; 


&process_cmd($run_Trinity_cmd, "$checkpoints_dir/trinity.ok");

# Examine top of Trinity.fasta file
&process_cmd("head trinity_out_dir/Trinity.fasta", "$checkpoints_dir/head_trinity.ok");

# count the number of transcripts assembled.
&process_cmd("grep '>' trinity_out_dir/Trinity.fasta | wc -l ", "$checkpoints_dir/count_trans.ok");


# Get Trinity stats:
&process_cmd("$trinity_dir/util/TrinityStats.pl trinity_out_dir/Trinity.fasta", "$checkpoints_dir/trin_stats.ok");


## representation of reads by the assembly
&process_cmd("$trinity_dir/util/bowtie_PE_separate_then_join.pl --target trinity_out_dir/Trinity.fasta --seqType fq --left data/wt_SRR1582651_1.fastq --right data/wt_SRR1582651_2.fastq --aligner bowtie -- -p 2 --all --best --strata -m 300", 
             "$checkpoints_dir/bowtie_PE_sep_join.ok");

&process_cmd("ls -ltr bowtie_out/", "$checkpoints_dir/ls_bowtie_outdir.ok"); 

&process_cmd("$trinity_dir/util/SAM_nameSorted_to_uniq_count_stats.pl bowtie_out/bowtie_out.nameSorted.bam",
             "$checkpoints_dir/nameSorted_sam_stats.ok");


###########################################
## assess number of full-length transcripts

&process_cmd("blastx -query trinity_out_dir/Trinity.fasta -db data/mini_sprot.pep -out blastx.outfmt6 -evalue 1e-20 -num_threads 2 -max_target_seqs 1 -outfmt 6", "$checkpoints_dir/blastx_for_full_length.ok");

&process_cmd("$trinity_dir/util/analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 trinity_out_dir/Trinity.fasta data/mini_sprot.pep", "$checkpoints_dir/tophat_blast_cov_stats.ok");


###################################
## Abundance estimation using RSEM
###################################

# also, write the samples.txt file

open (my $ofh, ">samples.txt") or die "Error, cannot write to file samples.txt";
my @rsem_trans_result_files;
my @rsem_gene_result_files;
my @bam_files;


my $first_round = 1;
foreach my $condition (sort keys %samples) {
    
    my $samples_href = $samples{$condition};
    
    foreach my $sample (keys %$samples_href) {
        
        print $ofh join("\t", $condition, $sample) . "\n";
        
        my $replicates_aref = $samples_href->{$sample};
        
        my ($left_fq, $right_fq) = @{$replicates_aref};
        
        my $output_dir = "$sample.RSEM";
        
        my $rsem_trans_result_file = "$output_dir/$sample.isoforms.results";
        push (@rsem_trans_result_files, $rsem_trans_result_file);
        
        my $rsem_gene_result_file = "$output_dir/$sample.genes.results";
        push (@rsem_gene_result_files, $rsem_gene_result_file);
        
        
        my $align_estimate_command = "$trinity_dir/util/align_and_estimate_abundance.pl --seqType fq "
            . " --left $left_fq --right $right_fq "
            . " --transcripts trinity_out_dir/Trinity.fasta "
            . " --output_prefix $sample --est_method RSEM "
            . " --aln_method bowtie --trinity_mode --prep_reference --coordsort_bam "
            . " --output_dir $output_dir";
        
        &process_cmd($align_estimate_command, "$checkpoints_dir/$sample.align_estimate.ok");
        
        push (@bam_files, "$output_dir/$sample.bowtie.csorted.bam");
        
        if ($first_round) {
            # look at the output
            &process_cmd("head $rsem_trans_result_file", "$checkpoints_dir/head.$sample.rsem.trans.ok");
            
            &process_cmd("head $rsem_gene_result_file", "$checkpoints_dir/head.$sample.rsem.genes.ok");
            
        }
    }
    
    
}
close $ofh; # samples.txt

## generate matrix of counts and perform TMM normalization
&process_cmd("$trinity_dir/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix Trinity_trans @rsem_trans_result_files", "$checkpoints_dir/trans_matrices.ok");

## Look at the counts matrix
&process_cmd("head -n20 Trinity_trans.counts.matrix", "$checkpoints_dir/head.counts.matrix.ok");

## Look at the expression matrix:
&process_cmd("head -n20 Trinity_trans.TMM.EXPR.matrix", "$checkpoints_dir/head.expr.matrix.ok");

# make the gene matrices
&process_cmd("$trinity_dir/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix Trinity_genes @rsem_gene_result_files", "$checkpoints_dir/gene_matrices.ok");

&process_cmd("ls -1 | grep gene | grep matrix", "$checkpoints_dir/ls_gene_matrices.ok");

## Examine the E90N50 statistic
&process_cmd("$trinity_dir//util/misc/contig_ExN50_statistic.pl  Trinity_trans.TMM.EXPR.matrix trinity_out_dir/Trinity.fasta > ExN50.stats", "$checkpoints_dir/ExNstats.ok");

&process_cmd("cat ExN50.stats", "$checkpoints_dir/cat_ExNstats.ok");

## Plot the values
&process_cmd("$trinity_dir/util/misc/plot_ExN50_statistic.Rscript ExN50.stats", "$checkpoints_dir/plot_ExN50.ok");

&show("ExN50.stats.plot.pdf");

## Examine read alignments in IGV

my $igv_cmd = "igv.sh -g trinity_out_dir/Trinity.fasta " . join(",", @bam_files);
if ($AUTO_MODE) {
    $igv_cmd .= " & ";
}
&process_cmd($igv_cmd, "$checkpoints_dir/igv_trinity_reads.ok");


##############
## DE analysis
##############

## make a samples file
&process_cmd("cat samples.txt", "$checkpoints_dir/examine_samples_txt.ok");

&process_cmd("cat -te samples.txt", "$checkpoints_dir/example_samples_txt_cat_te.ok");

## run edgeR
&process_cmd("$trinity_dir/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Trinity_trans.counts.matrix --samples_file samples.txt --method edgeR --output edgeR", "$checkpoints_dir/run.edgeR.ok");

# take a look at what edgeR generated:
&process_cmd("ls -ltr edgeR_trans/", "$checkpoints_dir/ls.edgeR_trans.dir.ok");


&process_cmd("head edgeR_trans/Trinity_trans.counts.matrix.GSNO_vs_WT.edgeR.DE_results", "$checkpoints_dir/head.edgeR_trans.DE_results.ok");

&show("edge_transR/Trinity_trans.counts.matrix.GSNO_vs_WT.edgeR.DE_results.MA_n_Volcano.pdf");

&change_dir("edgeR_trans", "$checkpoints_dir/cd.edgeR.ok");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../Trinity_trans.TMM.EXPR.matrix --samples ../samples.txt -P 1e-3 -C 2",
             "$checkpoints_dir/analyze_diff_expr.ok");


&process_cmd("wc -l diffExpr.P1e-3_C2.matrix", "$checkpoints_dir/wc_diff_expr_matrix.ok"); # number of DE transcripts + 1

&show("diffExpr.P1e-3_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData",
             "$checkpoints_dir/cut_clusters_tree.ok");

&show("diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60/my_cluster_plots.pdf");


## Now run the DE analysis for the genes

&change_dir("../", "$checkpoints_dir/cd_back_to_wd_after_edgeR.ok");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Trinity_genes.counts.matrix --samples_file samples.txt  --method edgeR --output edgeR_gene", "$checkpoints_dir/gene_DE_analysis.ok");

&process_cmd("ls -ltr edgeR_gene/", "$checkpoints_dir/ls_edgeR_gene_dir.ok");

#########################
## Time now for Trinotate
#########################

&run_cmd("pwd", "$checkpoints_dir/ensure_pwd_pre_Trinotate.ok");


&run_Trinotate_demo(); # cd's into Trinotate

################
## GO enrichment 
################

&change_dir("../edgeR_trans", "$checkpoints_dir/cd_back_to_edgeR_after_trinotate.ok");

## need sequence lengths file
&process_cmd("$trinity_dir/util/misc/fasta_seq_length.pl ../trinity_out_dir/Trinity.fasta > Trinity.seqLengths", "$checkpoints_dir/trin_seqlengths.ok");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/run_GOseq.pl --genes_single_factor Trinity_trans.counts.matrix.GSNO_vs_WT.edgeR.DE_results.P1e-3_C2.GSNO-UP.subset --GO_assignments ../Trinotate/Trinotate.xls.gene_ontology --lengths Trinity.seqLengths", "$checkpoints_dir/go_seq_gsno.ok");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/run_GOseq.pl --genes_single_factor Trinity_trans.counts.matrix.GSNO_vs_WT.edgeR.DE_results.P1e-3_C2.WT-UP.subset --GO_assignments ../Trinotate/Trinotate.xls.gene_ontology --lengths Trinity.seqLengths", "$checkpoints_dir/go_seq_wt.ok");


#######################################
## Prep TrinotateWeb w/ Expression Data
#######################################

&change_dir("../Trinotate", "$checkpoints_dir/cd_back_to_Trinotate_from_edgeR.ok");

&process_cmd("$trinotate_dir/util/transcript_expression/import_expression_and_DE_results.pl --sqlite Trinotate.sqlite --transcript_mode --samples_file ../samples.txt --count_matrix ../Trinity_trans.counts.matrix --fpkm_matrix ../Trinity_trans.TMM.EXPR.matrix",
             "$checkpoints_dir/Trinotate.load_trans_expr_data.ok");

&process_cmd("$trinotate_dir/util/transcript_expression/import_expression_and_DE_results.pl --sqlite Trinotate.sqlite --transcript_mode --samples_file ../samples.txt --DE_dir ../edgeR_trans",
    "$checkpoints_dir/Trinotate.load_trans_DE_data.ok");

&process_cmd("$trinotate_dir/util/transcript_expression/import_transcript_clusters.pl --sqlite Trinotate.sqlite --group_name DE_all_vs_all --analysis_name diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60 ../edgeR_trans/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60/*matrix",
             "$checkpoints_dir/Trinotate.load_expr_clusters.ok");

# load in the gene results

&process_cmd("$trinotate_dir/util/transcript_expression/import_expression_and_DE_results.pl  --sqlite ../Trinotate.sqlite --component_mode  --samples_file ../samples.txt --count_matrix ../Trinity_genes.counts.matrix --fpkm_matrix ../Trinity_genes.TMM.EXPR.matrix", "$checkpoints_dir/Trinotate.load_gene_expr_data.ok");

&process_cmd("$trinotate_dir/util/transcript_expression/import_expression_and_DE_results.pl --sqlite ../Trinotate.sqlite --component_mode  --samples_file ../samples.txt --DE_dir ../edgeR_gene ", "$checkpoints_dir/Trinotate.load_gene_DE_data.ok");


print STDERR "\n\n\tCommand-line Demo complete.  Congratulations! :)  Now explore your data via TrinotateWeb\n\n\n\n";


exit(0);

####
sub process_cmd {
    my ($cmd, $checkpoint) = @_;

    unless ($checkpoint) {
        die "Error, need checkpoint file defined";
    }
    
    $STEP_COUNT++;
    $checkpoint .= "s_" . sprintf("%02i", $STEP_COUNT);
    

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


####
sub change_dir {
    my ($dest_dir, $checkpoint) = @_;

    
    eval {
        &process_cmd("cd $dest_dir", $checkpoint);
    };
    if ($@) {
        print STDERR "\n ** Note: if you see an error message above about not being able to cd, just ignore it... it's a weird behavior of this demo script. Rest assured we've \'cd $dest_dir\' just fine.   :)\n\n";
        system("touch $checkpoint");
    }

    # now do it in the script. :)
    chdir("$dest_dir") or die "Error, could not cd to $dest_dir"; 
 
    return;
}

####
sub run_Trinotate_demo {
   
    &process_cmd("mkdir Trinotate", "$checkpoints_dir/mkdir_Trinotate.ok");

    &change_dir("Trinotate", "$checkpoints_dir/cd_Trinotate.ok");
    
    my $ini_reader = new IniReader($trinotate_conf_file);
    
    my @sections = $ini_reader->get_section_headings();
    @sections = grep { $_ ne 'GLOBALS' } @sections;

    my %globals = $ini_reader->get_section_hash('GLOBALS');
    $globals{TRANSCRIPTS_FASTA} = "../trinity_out_dir/Trinity.fasta";
    $globals{GENE_TO_TRANS_MAP} = "../trinity_out_dir/Trinity.fasta.gene_trans_map";
    $globals{CPU} = 2;
    $globals{TRINOTATE_HOME} = $trinotate_dir;
    
    ## get command structs
    my @cmd_structs;
    foreach my $section (@sections) {
        my %keyvals = $ini_reader->get_section_hash($section);
        $keyvals{__SECTION__} = $section;
        
        if ($keyvals{RUN} =~ /^T/i) {
            push (@cmd_structs, \%keyvals);
        }
    }

    @cmd_structs = sort {$a->{RANK}<=>$b->{RANK}} @cmd_structs;

    foreach my $cmd_struct (@cmd_structs) {
        my $CMD = $cmd_struct->{CMD};
        $CMD = &substitute_tokens($CMD, \%globals);
        
        my $section_name = $cmd_struct->{__SECTION__};
        my $checkpoint_file = "$checkpoints_dir/Trinotate.$section_name.ok";
        
        &process_cmd($CMD, $checkpoint_file);
        
    }
    
    return;
}

####
sub substitute_tokens {
    my ($cmd, $globals_href) = @_;

    my %token_templates;
    while ($cmd =~ /(\{__\S+__\})/g) {
        my $token_template = $1;
        
        $token_templates{$token_template}++;
    }

    if (%token_templates) {
        foreach my $token_template (keys %token_templates) {
            $token_template =~ /\{__(\S+)__\}/ or die "Error, not able to parse token template: $token_template";
            my $token_name = $1;

            my $replacement_val = $globals_href->{$token_name};
            unless (defined $replacement_val) {
                die "Error, unable to identify global value for token name: $token_name of cmd: $cmd";
            }
            $cmd =~ s/$token_template/$replacement_val/g;
        }
    }

    return($cmd);
}
