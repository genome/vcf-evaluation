#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

my $DEBUG = 0;
my $VCFLIB= "/gscmnt/gc3042/cle_validation/src/vcflib/bin";
my $SEED = '777';

# This script will look at the evaluation output directory
# For both the 0-blah and 1-blah
# It will extract file-unique sites
# It will subsample the subsequent VCF to approximately the correct variant number
# It will bgzip and tabix each VCF.
# This will give us VCFs ready for loading into IGV.


my $evaluation_directory;
my $variant_number;
my $output_directory;
my $help;

sub print_help {
    print STDERR __FILE__, " --evaluate-vcf-outdir <path> --variant-number <int> --outdir <path>\n";
    exit 1;
}

GetOptions(
    'evaluate-vcf-outdir=s' => \$evaluation_directory,
    'variant-number=i' => \$variant_number,
    'outdir=s' => \$output_directory,
    'help!' => \$help,
) or print_help();
print_help() if($help || !(defined $evaluation_directory && defined $variant_number && $output_directory));

my ($evaluated_vcf, $standard_vcf) = find_vcfs($evaluation_directory);

my $subsample_eval_vcf_output = "$output_directory/evaluated_only.vcf.gz";
my $subsample_std_vcf_output = "$output_directory/standard_only.vcf.gz";
subsample_vcf(':1,0:',$evaluated_vcf, $subsample_eval_vcf_output, $variant_number);
subsample_vcf(':0,1:',$standard_vcf, $subsample_std_vcf_output, $variant_number);

exit;

sub find_vcfs {
    my ($eval_dir) = @_;
    #note that these have gz appended but are NOT zipped
    my @evaluated_vcfs = glob($eval_dir . "/0-*vcf*");
    my @standard_vcfs = glob($eval_dir . "/1-*vcf*");
    if(scalar(@evaluated_vcfs) > 1) {
        die "Found more than one evaluated VCF file: " . join(" ", @evaluated_vcfs) . "\n";
    }
    unless(@evaluated_vcfs) {
        die "Failed to find evaluated VCF file\n";
    }
    if(scalar(@standard_vcfs) > 1) {
        die "Found more than one standard VCF file: " . join(" ", @standard_vcfs) . "\n";
    }
    unless(@standard_vcfs) {
        die "Failed to find standard VCF file\n";
    }
    #Above checks should implicitly ensure that we only have found one file.
    return (@evaluated_vcfs, @standard_vcfs);
}
    
sub subsample_vcf {
    my ($pattern, $vcf, $outname, $number) = @_;
    my $selection_regex_string = '^#\|' . $pattern;
    my $count = count_variants($pattern, $vcf);
    my $ratio = calculate_ratio($number, $count);
    #execute the selection
    #tabix index that bad boy
    my $cmd = qq{grep '$selection_regex_string' $vcf | $VCFLIB/vcfrandomsample -r $ratio -p $SEED | bgzip -c > $outname};
    execute($cmd);
    execute("tabix -p vcf $outname");
}

sub count_variants {
    my ($pattern, $vcf) = @_;
    my @results = `grep '$pattern' $vcf | grep -v '^#' | wc -l`;
    chomp $results[0];
    return $results[0];
}

sub calculate_ratio {
    my ($target, $count) = @_;
    if($count <= $target) {
        return 1;
    }
    else {
        return $target / $count;
    }
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
