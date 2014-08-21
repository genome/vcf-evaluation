#!/usr/bin/env perl

use Test::More;
use Test::Exception;
use File::Temp;

use strict;
use warnings;

use_ok("VcfCompare");

my $file_contents = <<'TESTFILE';
ours	theirs	type	some_sample
0	1	partial_match	369
1	0	partial_match	316
1	1	partial_match	16816
0	1	exact_match	375
1	0	exact_match	317
1	1	exact_match	16810
0	1	partial_miss	0
1	0	partial_miss	5
1	1	partial_miss	0
0	1	complete_miss	369
1	0	complete_miss	311
1	1	complete_miss	0
TESTFILE

my $tfh = File::Temp->new() or die "Unable to create tempfile\n";
print $tfh $file_contents;
$tfh->flush;


my $table = new_ok("VcfCompare", [$tfh->filename]);
is($table->unique_count("theirs", "partial_match", "some_sample"), 369, "test unique count for second file");
is($table->unique_count("ours", "partial_miss", "some_sample"), 5, "test unique count for first file");
dies_ok { $table->unique_count("yours", "partial_match", "some_sample") } "unique_count dies on invalid filename";
dies_ok { $table->unique_count("ours", "total_match", "some_sample") } "unique_count dies on invalid type";
dies_ok { $table->unique_count("ours", "partial_match", "some_samples") } "unique_count dies on invalid sample name";

is($table->joint_count("exact_match", "some_sample"), 16810, "test joint count");
is($table->joint_count("complete_miss", "some_sample"), 0, "test joint count for last line");
dies_ok { $table->unique_count("total_match", "some_sample") } "joint_count dies on invalid type";
dies_ok { $table->unique_count("partial_match", "some_samples") } "joint_count dies on invalid sample name";


done_testing();
