#!/usr/bin/perl

use strict;
use warnings;

use lib "/gscmnt/gc3042/cle_validation/src/genome/lib/perl/";

use File::Spec;
use Getopt::Long;
use Genome;
use Cwd;

my $model_group = shift;

for my $build (Genome::ModelGroup->get($model_group)->builds) {
    for my $variant_type qw( snvs indels) {

        my $file = $build->data_directory . "/variants/$variant_type.vcf.gz";
        my $sample;
        if($build->model->subclass_name eq 'Genome::Model::SomaticValidation') {
            $sample = $build->model->tumor_sample->name;
        }
        else {
            $sample = $build->model->subject_name;
        }
        print join("\t", $build->model->name, $build->id, $file, $variant_type, $sample), "\n";
    }
}
