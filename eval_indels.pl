#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome;

my $true_negative_bed = "/gscmnt/sata831/info/medseq/dlarson/cap_validation/mixing_experiment_with_platinum/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.minus_indels.targeted.sorted.bed";
my $true_positive_bed = "/gscmnt/sata831/info/medseq/dlarson/cap_validation/mixing_experiment_with_platinum/plat_indels_intersect_nist.on_target.bed";
my $target_region_bed = "/gscmnt/sata831/info/medseq/dlarson/cap_validation/mixing_experiment_with_platinum/VCRome.targets.sorted.bed";

my $som_var_model_group = shift;

my $true_negative_space = bed_size($true_negative_bed);

for my $build (Genome::ModelGroup->get($som_var_model_group)->builds) {
    my $dirname = $build->id;

    unless(-d $dirname) {
        mkdir $dirname;
    }

    my $on_target_snvs = $dirname . "/indels.hq.on_target.bed";
    my $fixed_up_indels = $dirname . "/indels.hq.on_target.uniform_null.bed";
    my $fn_bed = $dirname . "/tp_indel.missed";
    my $tp_bed = $dirname . "/tp_indel.found";
    #my $tn_bed = $dirname . "/tn_indel.not_called";
    my $fp_bed = $dirname . "/non_tp_indel.called";
    my $fp_in_tn_bed = $dirname . "/non_tp_indel.called.in_tn_bed";
    my $fp_outside_tn_bed = $dirname . "/non_tp_indel.called.outside_tn_bed";
    my $tp_bed_double_check = $dirname . "/tp_indel.reverse_intersect";

    generate_on_target_bed($build, $on_target_snvs);
    
    `sed 's/0\\//*\\//' $on_target_snvs | sed 's/\\/0/\\/*/' > $fixed_up_indels`;
    
    my ($true_positives, $false_negatives) = exact_compare($true_positive_bed, $fixed_up_indels, $fn_bed, $tp_bed);
    my ($true_positive_doublecheck, $false_positives_total_target ) = exact_compare($fixed_up_indels, $true_positive_bed, $fp_bed, $tp_bed_double_check);
    my ($false_positives, $false_positives_outside_true_negative) = compare($fp_bed, $true_negative_bed, $fp_outside_tn_bed, $fp_in_tn_bed);
    unless($true_positives == $true_positive_doublecheck) {
        die "True Positive variant number depended on order passed to Joinx: $true_positives != $true_positive_doublecheck\n";
    }

    print $build->model->name," (" . $build->id . ")","\t",join("\t", $true_positives, $false_negatives, $true_negative_space - $false_positives, $false_positives);
    print "\n";
}


sub snv_bed {
    my $build = shift;
    return $build->data_directory . "/variants/indels.hq.bed";
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

sub bed_size {
    my $bed = shift;
    my @results = `gmt bed sum-intervals --bed-file $bed`; #docs for sum-intervals incorrectly describe its operation. This WILL do the right thing currently.
    chomp $results[0];
    return $results[0];
}
