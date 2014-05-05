#!/usr/bin/perl

use warnings;
use strict;

use IO::File;
use Getopt::Long;
use File::Basename;
use File::Spec;

my $JOINX = "/usr/bin/joinx1.8";
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
my $help;

GetOptions(
    'vcf=s' => \$vcf,
    'roi=s' => \$roi,
    'gold-vcf=s' => \$gold_vcf,
    'true-negative-bed=s' => \$tn_bed,
    'help!' => \$help,
) or print_help();
print_help() if $help;
my $bgzip_pipe_cmd = "| bgzip -c ";

my ($basename, $path, $suffix) = fileparse($vcf, ".vcf.gz");

restrict("$basename$suffix", $roi, "$basename.roi.vcf.gz");
restrict($gold_vcf, $roi, "$gold_vcf.roi.vcf.gz");
restrict($tn_bed, $roi, "$tn_bed.roi.bed.gz");

pass_only("$basename.roi.vcf.gz", "$basename.roi.pass_only.vcf.gz");
allelic_primitives("$basename.roi.pass_only.vcf.gz", "$basename.roi.pass_only.allelic_primitives.vcf.gz");
normalize_vcf("$basename.roi.pass_only.allelic_primitives.vcf.gz", $REFERENCE, "$basename.roi.pass_only.allelic_primitives.normalized.vcf.gz");
sort_file("$basename.roi.pass_only.allelic_primitives.normalized.vcf.gz","$basename.roi.pass_only.allelic_primitives.normalized.sorted.vcf.gz");
restrict("$basename.roi.pass_only.allelic_primitives.normalized.sorted.vcf.gz", $roi, "$basename.roi.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz");
compare("$basename.roi.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz", "$gold_vcf.roi.vcf.gz", "$basename.roi.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz.compared");
number_within_roi("$basename.roi.pass_only.allelic_primitives.normalized.sorted.reroi.vcf.gz", "$tn_bed.roi.bed.gz", "$basename.roi.pass_only.allelic_primitives.normalized.sorted.reroi.in_tn_bed.vcf.gz");

sub print_help {
    print STDERR "evaluate_vcf --vcf --roi --gold-vcf --true-negative-bed\n";
    exit;
}

sub allelic_primitives {
    my ($input_file, $output_file) = @_;
    execute("$VCFLIB/vcfallelicprimitives -t ALLELICPRIMITIVE $input_file | $VCFLIB/vcffixup - |  $VCFLIB/vcffilter -f 'AC > 0' $bgzip_pipe_cmd > $output_file");
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

sub restrict {
    my ($input_file, $roi_file, $output_file)  = @_;

    #TODO Check on what happens with headers if $input_file has a header
    #TODO Check on what happens to VCF entries that span a boundary of the ROI (e.g. deletion)
    execute("zgrep '^#' $input_file | perl -pe 's/\\s(\\S+?NA12878)(.+?)(\\s+)/\\tNA12878\$3/g' > /tmp/header");
    my $cmd = "$BEDTOOLS intersect -a $input_file -b $roi_file | cat /tmp/header - $bgzip_pipe_cmd > $output_file";
    execute($cmd); #this is not very safe. I would really prefer to use Genome or IPC::Run
    execute("tabix -p vcf $output_file");
}

sub pass_only {
    my ($input_file, $output_file) = @_;

    execute("$VCFLIB/vcffilter -g 'FT = PASS | FT = .' $input_file | grep -v '.	.\$' $bgzip_pipe_cmd > $output_file");
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

sub true_positives {
    my ($joinx_output) = @_;
    my $shared_count = 0;
    my $gold_only_count = 0;
    my $test_only_count = 0;
    my $fh = IO::File->new($joinx_output) or die "Unable to open $joinx_output to calculate size\n";
    my $header = $fh->getline;
    my @lines = map { chomp $_; $_} $fh->getlines;
    $fh->close;

    $test_only_count += (split /\t/, $lines[0])[2];
    $gold_only_count += (split /\t/, $lines[1])[2];
    $shared_count += (split /\t/, $lines[2])[2];
    return ($test_only_count, $gold_only_count, $shared_count);
}


sub number_within_roi {
    my ($input_file, $roi, $output_file) = @_;
    execute("$BEDTOOLS intersect -header -a $input_file -b $roi $bgzip_pipe_cmd > $output_file");
    execute("tabix -p vcf $output_file");
    return count($output_file);
}


sub count {
    my $file = shift;
    my @results = `zgrep -v '^#' $file | cut -f1,2,3 | sort -u | wc -l`;
    chomp $results[0];
    return $results[0];
}

sub bed_size {
    my $bed = shift;
    my $count = 0;
    my $fh = IO::File->new($bed) or die "Unable to open $bed to calculate size\n";
    while(my $line  = $fh->getline) {
        chomp $line;
        my ($chr, $start, $stop) = split "\t", $line;
        $count += $stop-$start;
    }
    $fh->close;
    return $count;
}

sub print_versions {
    print `$JOINX --version`;
    print `$BEDTOOLS --version`;
}

sub execute {
    my $cmd = shift;

    print STDERR $cmd,"\n";
    if(!$DEBUG) {
        print `$cmd`;
    }
}
