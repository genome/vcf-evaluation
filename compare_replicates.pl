#!/usr/bin/perl

use warnings;
use strict;

use IO::File;
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Temp qw( tempfile );
use VcfCompare;

my $JOINX = "/gscmnt/gc3042/cle_validation/src/joinx/build/bin/joinx";
my $BEDTOOLS = "/gscmnt/gc3042/cle_validation/src/bedtools2-2.19.1/bin/bedtools";
my $VCFLIB= "/gscmnt/gc3042/cle_validation/src/vcflib/bin";
my $REFERENCE = "/gscmnt/gc3042/cle_validation/reference/all_sequences.fa";
my $DEBUG = 0;

#Basic algorithm should be as follows:
# 1. Take in the ROI to restrict to, the VCF to evaluate, a BED file of positions you shouldn't find and a VCF file of variants you SHOULD find. (note that your VCF shouldn't contain hom ref calls then)
# 2. For the VCF to be evaluated:
#   2.1. Restrict to ROI
#   2.2. Run allelic primitives
#   2.3. Normalize
#   2.4. Sort
#   2.5. Re-restrict to ROI
#   2.6. Run Joinx compare-gt to get the intersection. (I don't believe this supports partial overlap. That may need to be requested.)
#   2.7. Run bedtools intersect of the evaluation VCF against the TN BED file to get #FP
#   2.8. Count the input files
#   2.9. Report the numbers
my $vcf;
my $roi;
my $baseline_vcf;
my $old_sample;
my $new_sample;
my $baseline_sample;
my $clean_indels = 0;
my $pass_only_exp = q{-g 'FT = PASS | FT = .'};
my $help;

GetOptions(
    'vcf=s' => \$vcf,
    'old-sample=s' => \$old_sample,
    'new-sample=s' => \$new_sample,
    'roi=s' => \$roi,
    'baseline-vcf=s' => \$baseline_vcf,
    'baseline-sample=s' => \$baseline_sample,
    'pass-only-expression=s' => \$pass_only_exp,
    'clean-indels' => \$clean_indels,
    'help!' => \$help,
) or print_help();
print_help() if $help;
my $bgzip_pipe_cmd = "| bgzip -c ";

$vcf = clean_caf_and_bgzip($vcf);
$baseline_vcf = clean_caf_and_bgzip($baseline_vcf);

my $prepped_replicate_vcf = prepare_vcf($vcf, $roi, $clean_indels, "replicate", $pass_only_exp, $old_sample);
my $prepped_baseline_vcf = prepare_vcf($baseline_vcf, $roi, $clean_indels, "baseline", $pass_only_exp, $baseline_sample);
compare_partial($prepped_replicate_vcf, ".", $prepped_baseline_vcf, "replicate.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz.compared", $baseline_sample, $old_sample, $new_sample);

my %results = true_positives($prepped_replicate_vcf, $prepped_baseline_vcf, "replicate.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz.compared", $new_sample);

print join("\t", 
    $results{true_positive_exact}, 
    $results{true_positive_exact} + $results{false_negative_exact},
    ($results{true_positive_exact} + $results{false_negative_exact}) ? $results{true_positive_exact} / ($results{true_positive_exact} + $results{false_negative_exact}) : 'nan',
    $results{true_positive_partial},
    $results{true_positive_partial} + $results{false_negative_partial},
    ($results{true_positive_partial} + $results{false_negative_partial}) ? $results{true_positive_partial} / ($results{true_positive_partial} + $results{false_negative_partial}) : 'nan',
    $results{false_positive_exact},
    $results{false_positive_partial},
    ($results{false_positive_exact} + $results{true_positive_exact}) ? $results{true_positive_exact} / ($results{false_positive_exact} + $results{true_positive_exact}) : 'nan',
    ($results{false_positive_partial} + $results{true_positive_partial}) ? $results{true_positive_partial} / ($results{false_positive_partial} + $results{true_positive_partial}) : 'nan',
), "\n"; 


sub print_help {
    print STDERR "compare_replicates --vcf --roi --baseline-vcf \n";
    exit;
}

sub allelic_primitives {
    my ($input_file, $output_file) = @_;
    if(count($input_file)) {
        execute("$VCFLIB/vcfallelicprimitives -t ALLELICPRIMITIVE $input_file | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f 'AC > 0' $bgzip_pipe_cmd > $output_file");
    }
    else {
        execute("cp $input_file $output_file");
    }
    execute("tabix -p vcf $output_file");
}

sub sort_file {
    my ($input_file, $output_file) = @_;
    if(count($input_file)) {
        execute("$JOINX sort $input_file $bgzip_pipe_cmd > $output_file");
    }
    else {
        execute("cp $input_file $output_file");
    }
    execute("tabix -p vcf $output_file");
}

sub normalize_vcf {
    my ($input_file, $reference, $output_file) = @_;
    if(count($input_file)) {
        execute("$JOINX vcf-normalize-indels -f $reference $input_file $bgzip_pipe_cmd > $output_file");
    }
    else {
        execute("cp $input_file $output_file");
    }
    execute("tabix -p vcf $output_file");
}

sub restrict_vcf_to_sample {
    my ($input_file, $sample, $output_file) = @_;
    if(count($input_file)) {
        execute("$VCFLIB/vcfkeepsamples $input_file $sample $bgzip_pipe_cmd > $output_file");
    }
    else {
        execute("cp $input_file $output_file");
    }
    execute("tabix -p vcf $output_file");
}

sub restrict {
    my ($input_file, $roi_file, $output_file, $clean_indels)  = @_;

    #TODO Check on what happens with headers if $input_file has a header
    #TODO Check on what happens to VCF entries that span a boundary of the ROI (e.g. deletion)
    #NOTE bedtools intersect and vcflib vcfintersect both only look at the reference coordinate.
    #this means that deletions where the variant is outside of the ROI are pulled in
    #we don't want to do this so we will clean it up.
    if(count($input_file)) {
        if($clean_indels) {
            my ($tfh, $tempfilename) = tempfile();
            my $cmd = "zcat $input_file | $BEDTOOLS intersect -header -a stdin -b $roi_file > $tempfilename";
            execute($cmd); #this is not very safe. I would really prefer to use Genome or IPC::Run
            my %bed_ends;
            my $fh = IO::File->new($roi_file) or die "Unable to open BED file $roi_file for removal of bad indel lines\n";
            while(my $bedline = $fh->getline) {
                chomp $bedline;
                my ($chr, $start, $stop) = split "\t", $bedline;
                $bed_ends{"$chr\t$stop"} = 1;
            }
            $fh->close;
            $tfh->flush;
            my $ofh = IO::File->new("| bgzip -c > $output_file") or die "Unable to open output bgzip pipe\n";
            while(my $vcfline = $tfh->getline) {
                my ($chr, $pos) = split "\t", $vcfline;
                if($vcfline =~ /^#/ || !exists($bed_ends{"$chr\t$pos"})) {
                    print $ofh $vcfline;
                }
            }
            $ofh->close;
            $tfh->close;
        }
        else {
            my $cmd = "zcat $input_file | $BEDTOOLS intersect -header -a stdin -b $roi_file $bgzip_pipe_cmd > $output_file";
            execute($cmd); #this is not very safe. I would really prefer to use Genome or IPC::Run
        }
    }
    else {
        execute("cp $input_file $output_file");
    }
    execute("tabix -p vcf $output_file");
}

sub pass_only {
    my ($input_file, $output_file, $expression) = @_;
    if(count($input_file)) {

        if($expression) {
            execute("$VCFLIB/vcffilter $expression $input_file $bgzip_pipe_cmd > $output_file");
        }
        else {
            execute("zcat $input_file | perl -ape '\$_=q{} unless(\$F[6] eq q{PASS} || \$F[6] eq q{.} || \$F[0] =~ /^#/)' $bgzip_pipe_cmd > $output_file");
        }
    }
    else {
        execute("cp $input_file $output_file");
    }
    execute("tabix -p vcf $output_file");
}

sub compare {
#$ joinx1.8 vcf-compare-gt -h
#  -h [ --help ]            this message
#  -i [ --input-file ] arg  input file(s) (positional arguments work also)
#  -n [ --name ] arg        meaningful names for each of the input files (given
#                           in the same order)
#  -s [ --sample-name ] arg operate only on these samples (may specify multiple
#                           times)
    my ($input_file, $gold_file, $output_file) = @_;
    execute("$JOINX vcf-compare-gt $input_file $gold_file -s NA12878 > $output_file");
}

sub compare_partial {
    my ($input_file, $variant_directory, $gold_file, $output_file, $gold_sample, $eval_sample, $new_sample) = @_;
    my $rename_option = "";
    if($new_sample) {
        if($gold_sample) {
            $rename_option .= " -R $gold_sample=$new_sample";
        }
        if($eval_sample) {
            $rename_option .= " -R $eval_sample=$new_sample";
        }
    }
    execute("$JOINX vcf-compare $rename_option -d $variant_directory $input_file $gold_file -s $new_sample > $output_file");
}

sub true_positives {
    my ($input_file, $gold_file, $joinx_output, $new_sample) = @_;
    my $table = VcfCompare->new($joinx_output);
    #for now only do perfect matches
    return (
        false_positive_exact => $table->unique_count($input_file, "exact_match", $new_sample),
        false_negative_exact => $table->unique_count($gold_file, "exact_match", $new_sample),
        true_positive_exact => $table->joint_count("exact_match", $new_sample),
        false_positive_partial => $table->unique_count($input_file, "partial_match", $new_sample),
        false_negative_partial => $table->unique_count($gold_file, "partial_match", $new_sample),
        true_positive_partial => $table->joint_count("partial_match", $new_sample),
        false_positive_partial_miss => $table->unique_count($input_file, "partial_miss", $new_sample),
        false_negative_partial_miss => $table->unique_count($gold_file, "partial_miss", $new_sample),
        false_positive_complete_miss => $table->unique_count($input_file, "complete_miss", $new_sample),
        false_negative_complete_miss => $table->unique_count($gold_file, "complete_miss", $new_sample),
    );
}


sub number_within_roi {
    my ($input_file, $roi, $output_file) = @_;
    if(count($input_file)) {
        execute("zcat $input_file | $BEDTOOLS intersect -header -a stdin -b $roi $bgzip_pipe_cmd > $output_file");
    }
    else {
        execute("cp $input_file $output_file");
    }
    execute("tabix -p vcf $output_file");
    return count($output_file);
}

sub bed_size {
    my $bed = shift;
    my $count = 0;
    my $fh = IO::File->new("zcat $bed |") or die "Unable to open $bed to calculate size\n";
    while(my $line  = $fh->getline) {
        chomp $line;
        my ($chr, $start, $stop) = split "\t", $line;
        $count += $stop-$start;
    }
    $fh->close;
    return $count;
}

sub count {
    my $file = shift;
    my @results = `zgrep -v '^#' $file | cut -f1,2,3 | sort -u | wc -l`;
    chomp $results[0];
    return $results[0];
}


sub print_versions {
    print `$JOINX --version`;
    print `$BEDTOOLS --version`;
}

sub execute {
    my $cmd = shift;

    print STDERR $cmd,"\n";
    if(!$DEBUG) {
        my $rvalue = `$cmd`;
        unless(defined $rvalue) {
            die "Error running $cmd\n";
        }
        else {
            print $rvalue;
        }
    }
}

sub clean_caf_and_bgzip {
    my ($vcf) = @_;
    if($vcf =~ /\.vcf$/) {
        execute("perl -pe 's/CAF=.*?[;      ]//' $vcf $bgzip_pipe_cmd > $vcf.gz");
        return $vcf . ".gz";
    }
    else {
        return $vcf;
    }
}

sub prepare_vcf {
    my ($vcf, $roi, $clean_indels, $new_basename, $pass_only_exp, $original_sample) = @_;
    #should perform all the necessary ROI restriction, normalization, allelic primitives and renaming necessary.

    #this nonsense below should be moved into a function
    restrict($vcf, $roi, "$new_basename.roi.vcf.gz");
    restrict_vcf_to_sample("$new_basename.roi.vcf.gz", $original_sample, "$new_basename.roi.$original_sample.vcf.gz");

    pass_only("$new_basename.roi.$original_sample.vcf.gz", "$new_basename.roi.$original_sample.pass_only.vcf.gz", $pass_only_exp);
    allelic_primitives("$new_basename.roi.$original_sample.pass_only.vcf.gz", "$new_basename.roi.$original_sample.pass_only.allelic_primitives.vcf.gz");
    normalize_vcf("$new_basename.roi.$original_sample.pass_only.allelic_primitives.vcf.gz", $REFERENCE, "$new_basename.roi.$original_sample.pass_only.allelic_primitives.normalized.vcf.gz");
    sort_file("$new_basename.roi.$original_sample.pass_only.allelic_primitives.normalized.vcf.gz","$new_basename.roi.$original_sample.pass_only.allelic_primitives.normalized.sorted.vcf.gz");
    restrict("$new_basename.roi.$original_sample.pass_only.allelic_primitives.normalized.sorted.vcf.gz", $roi, "$new_basename.roi.$original_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz", $clean_indels); #we will clean indels here in case allelic primitives and indel normalization have moved things to the edges of ROIs and need to be cleaned
    return "$new_basename.roi.$original_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz";
}
