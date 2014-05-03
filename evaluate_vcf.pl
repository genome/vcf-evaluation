#!/usr/bin/perl

use warnings;
use strict;

use IO::File;

my $JOINX = "/usr/bin/joinx1.8";
my $BEDTOOLS = "/gscmnt/gc3042/cle_validation/src/bedtools2-2.19.1/bin/bedtools";
my $VCFLIB= "/gscmnt/gc3042/cle_validation/src/vcflib/bin";
my $REFERENCE = "/gscmnt/gc3042/cle_validation/reference/all_sequences.fa";

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
sub allelic_primitives {
    my ($input_file, $output_file) = @_;
    execute("$VCFLIB/vcfallelicprimitives -t ALLELICPRIMITIVE $input_file > $output_file");
}


sub sort {
    my ($input_file, $output_file) = @_;
    execute("$JOINX sort $input_file > $output_file");
}

sub normalize_vcf {
    my ($input_file, $reference, $output_file) = @_;
    execute("$JOINX vcf-normalize-indels -f $reference $input_file > $output_file");
}

sub restrict {
    my ($input_file, $roi_file, $output_file)  = @_;

    #TODO Check on what happens with headers if $input_file has a header
    #TODO Check on what happens to VCF entries that span a boundary of the ROI (e.g. deletion)
    my $cmd = "$BEDTOOLS subtract -a $input_file -b $roi_file > $output_file";
    execute($cmd); #this is not very safe. I would really prefer to use Genome or IPC::Run
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
    execute("$JOINX vcf-compare-gt $input_file $gold_file > $output_file");
}

sub number_within_roi {
    my ($input_file, $roi, $output_file) = @_;
    execute("$BEDTOOLS intersect -a -b > $output_file");
    count_bed($output_file);
}


sub count_bed {
    my $bed = shift;
    my @results = `cut -f1,2,3 $bed | sort -u | wc -l`;
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
