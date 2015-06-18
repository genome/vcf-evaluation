#!/usr/bin/perl

use strict;
use warnings;

use lib "/gscmnt/gc3042/cle_validation/src/genome/lib/perl/";

use File::Spec;
use File::Basename;
use Getopt::Long;
use Cwd;
use IO::File;

my $EVAL_SCRIPT_LOCATION = "/gscmnt/gc3042/cle_validation/src/evaluation/";
#inputs: reference, tn bed, indel variant vcf, snv variant vcf, model group, roi
my $config;
my $roi;
my $baseline_snv_vcf;
my $baseline_indel_vcf;
my $help;
my $baseline_sample;
my $expression;

GetOptions(
    'config=s' => \$config,
    'roi=s' => \$roi,
    'baseline-snv-vcf=s' => \$baseline_snv_vcf,
    'baseline-indel-vcf=s' => \$baseline_indel_vcf,
    'baseline-sample=s' => \$baseline_sample,
    'pass-filter-expression=s' => \$expression,
    'help!' => \$help,
) or print_help();
print_help() if $help;

$roi = Cwd::abs_path($roi);
$baseline_indel_vcf = Cwd::abs_path($baseline_indel_vcf);
$baseline_snv_vcf = Cwd::abs_path($baseline_snv_vcf);
$config = Cwd::abs_path($config);

print join("\t", qw( Name Id VarType True_Positive_Found_Exact Total_True_Positive_Exact Sensitivity_Exact  True_Positive_Found_Partial Total_True_Positive_Partial Sensitivity_Partial False_Positive_Exact False_Positive_Partial Exact_PPV Partial_PPV )),"\n";

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
    my $file_name = "replicate" . file_suffix($file);
    symlink $file, $file_name;
    my $roi_name = filename($roi);
    symlink $roi, $roi_name;
    my $baseline_file = $variant_type eq "snvs" ? $baseline_snv_vcf : $baseline_indel_vcf;
    my $baseline_file_name = "baseline" . file_suffix($baseline_file);
    symlink $baseline_file, $baseline_file_name;

    my $cmd = "perl -I $EVAL_SCRIPT_LOCATION $EVAL_SCRIPT_LOCATION/compare_replicates.pl --vcf $file_name --roi $roi_name --baseline-vcf $baseline_file_name --old-sample $sample --new-sample SAMPLE";
    if($baseline_sample) {
        $cmd .= " --baseline-sample $baseline_sample";
    }
    if(defined $expression) {
        $cmd .= qq{ --pass-only-expression "$expression"};
    }
    if($variant_type eq 'indels') {
        $cmd .= " --clean-indels";
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

sub file_suffix {
    my $path = shift;
    my ($basename, $fullpath, $suffix) = fileparse($path, ".vcf", ".vcf.gz",);
    return $suffix;
}

sub validate_type {
    my $type = shift;
    if($type ne 'snvs' && $type ne 'indels') {
        die "Invalid type: $type. Only snvs and indels are supported as valid types.\n";
    }
    return 1;   #return true
}
