#!/usr/bin/perl

use strict;
use warnings;

use lib "/gscmnt/gc3042/cle_validation/src/genome/lib/perl/";

use File::Spec;
use Getopt::Long;
use Genome;
use Cwd;

#inputs: reference, tn bed, indel variant vcf, snv variant vcf, model group, roi
my $model_group;
my $roi;
my $gold_snv_vcf;
my $gold_indel_vcf;
my $tn_bed;
my $help;

GetOptions(
    'model-group=s' => \$model_group,
    'roi=s' => \$roi,
    'gold-snv-vcf=s' => \$gold_snv_vcf,
    'gold-indel-vcf=s' => \$gold_indel_vcf,
    'true-negative-bed=s' => \$tn_bed,
    'help!' => \$help,
) or print_help();
print_help() if $help;

$roi = Cwd::abs_path($roi);
$gold_indel_vcf = Cwd::abs_path($gold_indel_vcf);
$gold_snv_vcf = Cwd::abs_path($gold_snv_vcf);
$tn_bed = Cwd::abs_path($tn_bed);

my $tn_bed_size;

for my $build (Genome::ModelGroup->get($model_group)->builds) {
    
    my $dirname = $build->id;
    unless(-d $dirname) {
        mkdir $dirname;
    }

    for my $variant_type qw( snvs indels ) {
        my $cwd = cwd();
        my $dirname = $build->id . "/$variant_type";
        unless(-d $dirname) {
            mkdir $dirname;
        }
        chdir $dirname;
        my $file = $build->data_directory . "/variants/$variant_type.vcf.gz";
        my $file_name = filename($file);
        symlink $file, $file_name;
        my $roi_name = filename($roi);
        symlink $roi, $roi_name;
        my $tn_bed_name = filename($tn_bed);
        symlink $tn_bed, $tn_bed_name;
        my $gold_file = $variant_type eq "snvs" ? $gold_snv_vcf : $gold_indel_vcf;
        my $gold_file_name = filename($gold_file);
        symlink $gold_file, $gold_file_name;

        my $old_sample = $build->model->subject_name;

        my @output = `perl /gscmnt/gc3042/cle_validation/src/evaluation/evaluate_vcf.pl --vcf $file_name --roi $roi_name --gold-vcf $gold_file_name --true-negative-bed $tn_bed_name --old-sample $old_sample --new-sample NA12878`;
        $tn_bed_size = bed_size($tn_bed_name . ".roi.bed.gz") unless defined $tn_bed_size;
        chomp($output[0]);
        print $build->model->name, "\t", $output[0], "\t", $tn_bed_size, "\n";
        chdir $cwd;
    }

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
