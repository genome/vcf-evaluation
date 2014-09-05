#!/usr/bin/perl

use strict;
use warnings;

use lib "/gscmnt/gc3042/cle_validation/src/genome/lib/perl/";

use File::Spec;
use Getopt::Long;
use Cwd;
use IO::File;

#inputs: reference, tn bed, indel variant vcf, snv variant vcf, model group, roi
my $config;
my $roi;
my $gold_snv_vcf;
my $gold_indel_vcf;
my $tn_bed;
my $help;
my $gold_sample;

GetOptions(
    'config=s' => \$config,
    'roi=s' => \$roi,
    'gold-snv-vcf=s' => \$gold_snv_vcf,
    'gold-indel-vcf=s' => \$gold_indel_vcf,
    'true-negative-bed=s' => \$tn_bed,
    'gold-sample=s' => \$gold_sample,
    'help!' => \$help,
) or print_help();
print_help() if $help;

$roi = Cwd::abs_path($roi);
$gold_indel_vcf = Cwd::abs_path($gold_indel_vcf);
$gold_snv_vcf = Cwd::abs_path($gold_snv_vcf);
$tn_bed = Cwd::abs_path($tn_bed);
$config = Cwd::abs_path($config);

my $tn_bed_size;
print join("\t", qw( Name Id VarType True_Positive_Found_Exact Total_True_Positive_Exact Sensitivity_Exact  True_Positive_Found_Partial Total_True_Positive_Partial Sensitivity_Partial False_Positive_Exact False_Positive_Partial True_Negatives Exact_Specificity Partial_Specificity Exact_PPV Partial_PPV VCF_Lines_Overlapping_TN Lines_Specificity_in_TN_Only)),"\n";

#will now read in a config file
#name\tid\tpath\ttype\tsample\n

my $fh = IO::File->new($config,"r") or die "Unable to open $config\n";

while(my $line = $fh->getline) {
    chomp $line;
    my ($name, $id, $path, $variant_type, $sample) = split "\t", $line;
    my $abs_path = Cwd::abs_path($path);
    validate_type($variant_type); #this will die if the type is not valid

    my $dirname = "${name}_${id}";
    unless(-d $dirname) {
        mkdir $dirname;
    }

    my $cwd = cwd();
    my $output_dirname = $dirname . "/$variant_type";
    unless(-d $output_dirname) {
        mkdir $output_dirname;
    }
    chdir $output_dirname;
    my $file = $abs_path;
    my $file_name = filename($file);
    symlink $file, $file_name;
    my $roi_name = filename($roi);
    symlink $roi, $roi_name;
    my $tn_bed_name = filename($tn_bed);
    symlink $tn_bed, $tn_bed_name;
    my $gold_file = $variant_type eq "snvs" ? $gold_snv_vcf : $gold_indel_vcf;
    my $gold_file_name = filename($gold_file);
    symlink $gold_file, $gold_file_name;

    #$tn_bed_size = bed_size($tn_bed_name . ".roi.bed.gz") unless defined $tn_bed_size;

    my $cmd = "perl -I ~dlarson/src/evaluate/ ~dlarson/src/evaluate/evaluate_vcf.pl --vcf $file_name --roi $roi_name --gold-vcf $gold_file_name --true-negative-bed $tn_bed_name --old-sample $sample --new-sample GOLDSTANDARD_SAMPLE";#--true-negative-size $tn_bed_size";
    if($gold_sample) {
        $cmd .= " --gold-sample $gold_sample";
    }
    print STDERR $cmd,"\n";
    my @output = `$cmd`;
    chomp($output[0]);
    print join("\t",$name, $id, $variant_type, $output[0]), "\n";
    chdir $cwd;
}


sub filename {
    my $path = shift;
    my @stuff = File::Spec->splitpath( $path );
    return $stuff[2];
}

sub print_help {
    print STDERR "evaluate_modelgroup --model-group --roi --gold-snv-vcf --gold-indel-vcf --true-negative-bed\n";
    exit;
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

sub validate_type {
    my $type = shift;
    if($type ne 'snvs' && $type ne 'indels') {
        die "Invalid type: $type. Only snvs and indels are supported as valid types.\n";
    }
    return 1;   #return true
}
