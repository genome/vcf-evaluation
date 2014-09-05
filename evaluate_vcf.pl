#!/usr/bin/perl

use warnings;
use strict;

use IO::File;
use Getopt::Long;
use File::Basename;
use File::Spec;
use VcfCompare;

my $JOINX = "~dlarson/src/joinx/build/bin/joinx";
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
my $gold_vcf;
my $tn_bed;
my $old_sample;
my $new_sample;
my $tn_bed_size;
my $gold_sample;
my $pass_only_exp = q{-g 'FT = PASS | FT = .'};
my $help;

GetOptions(
    'vcf=s' => \$vcf,
    'old-sample=s' => \$old_sample,
    'new-sample=s' => \$new_sample,
    'roi=s' => \$roi,
    'gold-vcf=s' => \$gold_vcf,
    'gold-sample=s' => \$gold_sample,
    'true-negative-bed=s' => \$tn_bed,
    'true-negative-size=i' => \$tn_bed_size,
    'pass-only-expression=s' => \$pass_only_exp,
    'help!' => \$help,
) or print_help();
print_help() if $help;
my $bgzip_pipe_cmd = "| bgzip -c ";

my ($basename, $path, $suffix) = fileparse($vcf, ".vcf.gz");

restrict("$basename$suffix", $roi, "$basename.roi.vcf.gz");
restrict($gold_vcf, $roi, "$gold_vcf.roi.vcf.gz");
restrict($tn_bed, $roi, "$tn_bed.roi.bed.gz");

restrict_vcf_to_sample("$basename.roi.vcf.gz", $old_sample, "$basename.roi.$old_sample.vcf.gz");

pass_only("$basename.roi.$old_sample.vcf.gz", "$basename.roi.$old_sample.pass_only.vcf.gz", $pass_only_exp);
allelic_primitives("$basename.roi.$old_sample.pass_only.vcf.gz", "$basename.roi.$old_sample.pass_only.allelic_primitives.vcf.gz");
normalize_vcf("$basename.roi.$old_sample.pass_only.allelic_primitives.vcf.gz", $REFERENCE, "$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.vcf.gz");
sort_file("$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.vcf.gz","$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.vcf.gz");
restrict("$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.vcf.gz", $roi, "$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz");
compare_partial("$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz", ".", "$gold_vcf.roi.vcf.gz", "$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz.compared", $gold_sample, $old_sample, $new_sample);

#NOTE We will not calculate the size of the roi here and instead will assume it is calculated elsewhere if needed.
$tn_bed_size = bed_size("$tn_bed.roi.bed.gz") unless defined $tn_bed_size;

my $false_positives_in_roi = number_within_roi("$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz", "$tn_bed.roi.bed.gz", "$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.reroi.in_tn_bed.vcf.gz");
my %results = true_positives("$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz", "$gold_vcf.roi.vcf.gz", "$basename.roi.$old_sample.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz.compared", $new_sample);
print join("\t", 
    $results{true_positive_exact}, 
    $results{true_positive_exact} + $results{false_negative_exact},
    $results{true_positive_exact} / ($results{true_positive_exact} + $results{false_negative_exact}),
    $results{true_positive_partial},
    $results{true_positive_partial} + $results{false_negative_partial},
    $results{true_positive_partial} / ($results{true_positive_partial} + $results{false_negative_partial}),
    $results{false_positive_exact},
    $results{false_positive_partial},
    $tn_bed_size,
    #exact specificity, these are not strictly accurate as the tn_bed may be significantly smaller than the target space ROI
    ($tn_bed_size - $results{false_positive_exact}) / $tn_bed_size,
    ($tn_bed_size - $results{false_positive_partial}) / $tn_bed_size,
    $results{true_positive_exact} / ($results{false_positive_exact} + $results{true_positive_exact}),
    $results{true_positive_partial} / ($results{false_positive_partial} + $results{true_positive_partial}),
    $false_positives_in_roi,
    ($tn_bed_size - $false_positives_in_roi) / $tn_bed_size, #this is more accurate than the other specificity measures, but doesn't take partial into account
), "\n"; 


sub print_help {
    print STDERR "evaluate_vcf --vcf --roi --gold-vcf --true-negative-bed\noptionally add --true-negative-size to avoid calculating the size of the true negative bed file by specifying it directly.\n";
    exit;
}

sub allelic_primitives {
    my ($input_file, $output_file) = @_;
    execute("$VCFLIB/vcfallelicprimitives -t ALLELICPRIMITIVE $input_file | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f 'AC > 0' $bgzip_pipe_cmd > $output_file");
    execute("tabix -p vcf $output_file");
}

sub sort_file {
    my ($input_file, $output_file) = @_;
    execute("$JOINX sort $input_file $bgzip_pipe_cmd > $output_file");
    execute("tabix -p vcf $output_file");
}

sub normalize_vcf {
    my ($input_file, $reference, $output_file) = @_;
    execute("$JOINX vcf-normalize-indels -f $reference $input_file $bgzip_pipe_cmd > $output_file");
    execute("tabix -p vcf $output_file");
}

sub restrict_vcf_to_sample {
    my ($input_file, $sample, $output_file) = @_;
    execute("$VCFLIB/vcfkeepsamples $input_file $sample $bgzip_pipe_cmd > $output_file");
    execute("tabix -p vcf $output_file");
}


sub restrict {
    my ($input_file, $roi_file, $output_file)  = @_;

    #TODO Check on what happens with headers if $input_file has a header
    #TODO Check on what happens to VCF entries that span a boundary of the ROI (e.g. deletion)
    my $replace_cmd = "";
    #if($old_sample && $new_sample) {
    #    $replace_cmd = "| perl -pe 's/$old_sample/$new_sample/g'";
    #}
    execute("zgrep '^#' $input_file $replace_cmd > /tmp/header");
    my $cmd = "zcat $input_file | $BEDTOOLS intersect -a stdin -b $roi_file | cat /tmp/header - $bgzip_pipe_cmd > $output_file";
    execute($cmd); #this is not very safe. I would really prefer to use Genome or IPC::Run
    execute("tabix -p vcf $output_file");
}

sub pass_only {
    my ($input_file, $output_file, $expression) = @_;

    if($expression) {
        execute("$VCFLIB/vcffilter $expression $input_file $bgzip_pipe_cmd > $output_file");
    }
    else {
        execute("zcat $input_file | perl -ape '\$_=q{} unless(\$F[6] eq q{PASS} || \$F[6] eq q{.} || \$F[0] =~ /^#/)' $bgzip_pipe_cmd > $output_file");
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
    execute("zcat $input_file | $BEDTOOLS intersect -header -a stdin -b $roi $bgzip_pipe_cmd > $output_file");
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
