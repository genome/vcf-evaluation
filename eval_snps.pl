#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome;

my $true_negative_bed = "/gscmnt/sata831/info/medseq/dlarson/cap_validation/mixing_experiment_with_platinum/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.minus_snps.targeted.1bp_window.no_variant.sorted.bed";
my $true_positive_bed = "/gscmnt/sata831/info/medseq/dlarson/cap_validation/mixing_experiment_with_platinum/plat_snps_intersect_nist.on_target.bed";
my $target_region_bed = "/gscmnt/sata831/info/medseq/dlarson/cap_validation/mixing_experiment_with_platinum/VCRome.targets.sorted.bed";

my $som_var_model_group = shift;

for my $build (Genome::ModelGroup->get($som_var_model_group)->builds) {
    my $dirname = $build->id;

    unless(-d $dirname) {
        mkdir $dirname;
    }

    my $on_target_snvs = $dirname . "/snvs.hq.on_target.bed";
    my $fn_bed = $dirname . "/tp.missed";
    my $tp_bed = $dirname . "/tp.found";
    my $tn_bed = $dirname . "/tn.not_called";
    my $fp_bed = $dirname . "/tn.called";

    generate_on_target_bed($build, $on_target_snvs);
    
    my ($true_positives, $false_negatives) = exact_compare($true_positive_bed, $on_target_snvs, $fn_bed, $tp_bed);
    my ($false_positives, $true_negatives) = compare($true_negative_bed, $on_target_snvs, $tn_bed, $fp_bed);

    print $build->model->name," (" . $build->id . ")","\t",join("\t", $true_positives, $false_negatives, $true_negatives, $false_positives);
    print "\n";
}


sub snv_bed {
    my $build = shift;
    return $build->data_directory . "/variants/snvs.hq.bed";
}

sub generate_on_target_bed {
    my $build = shift;
    my $output_file_name = shift;

    my $rv = Genome::Model::Tools::Joinx::Intersect->execute(
        input_file_a => snv_bed($build),
        input_file_b => $target_region_bed,
        output_file => $output_file_name,
        use_version => 1.8,
    );
    unless($rv) {
        die "Unable to limit $build->id snv bed to exome capture target region\n";
    }
}

sub exact_compare {
    my ($standard, $file, $missed_from_standard, $found_from_standard) = @_;
    #joinx1.8 intersect --exact-pos --exact-allele --iub-match  -a ../1kg_nist_intersect.on_target.bed -b snvs.hq.on_target.bed | cut -f1,2,3 | sort | uniq | wc -l
    my $rv = Genome::Model::Tools::Joinx::Intersect->execute(
        input_file_a => $standard,
        input_file_b => $file,
        exact_pos => 1,
        exact_allele => 1,
        iub_match => 1,
        miss_a_file => $missed_from_standard,
        output_file => $found_from_standard,
        use_version => 1.8
    );
    unless($rv) {
        die "unable to compare $standard to $file\n";
    }

    return (count_bed($found_from_standard), count_bed($missed_from_standard));
}

sub compare {
    my ($standard, $file, $missed_from_standard, $found_from_standard) = @_;
    # exact pos doesn't seem to be functioning properly (reporting no overlaps when there are clearly overlaps wtf) if exact pos is on.
    # without exact pos, we get the same results as bedtools intersect. Think this is ok.
    my $rv = Genome::Model::Tools::Joinx::Intersect->execute(
        input_file_a => $standard,
        input_file_b => $file,
        miss_a_file => $missed_from_standard,
        output_file => $found_from_standard,
        use_version => 1.8
    );
    unless($rv) {
        die "unable to compare $standard to $file\n";
    }

    return (count_bed($found_from_standard), count_bed($missed_from_standard));
}

sub count_bed {
    my $bed = shift;
    my @results = `cut -f1,2,3 $bed | sort -u | wc -l`;
    chomp $results[0];
    return $results[0];
}


