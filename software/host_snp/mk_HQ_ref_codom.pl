#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Convert ref_tag FASTA into the HQ_ref_codom table expected by
# codom_calling.pl. This converter is enzyme-agnostic: it preserves whatever
# tag length was generated upstream by Extract_cut_site.pl.

my $input  = 'ref_tag';
my $output = 'HQ_ref_codom';
my $help   = 0;
GetOptions(
    'i=s'    => \$input,
    'o=s'    => \$output,
    'h|help' => \$help,
) or usage(1);
usage(0) if $help;

open(my $IN,  '<', $input)  or die "Cannot read $input: $!\n";
open(my $OUT, '>', $output) or die "Cannot write $output: $!\n";

my ($id, $seq) = (undef, '');
my ($count, $expected_len) = (0, undef);
while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    if ($line =~ /^>/) {
        write_record($OUT, $id, $seq, \$count, \$expected_len) if defined $id;
        $id  = $line;
        $seq = '';
    } else {
        die "Sequence before FASTA header in $input\n" unless defined $id;
        $line =~ s/\s+//g;
        $seq .= uc($line);
    }
}
write_record($OUT, $id, $seq, \$count, \$expected_len) if defined $id;

close $IN;
close $OUT;
die "No reference tags were found in $input\n" unless $count;
print STDERR "Generated: $output ($count reference tags; tag length $expected_len bp)\n";

sub write_record {
    my ($OUT, $id, $seq, $count_ref, $len_ref) = @_;
    die "Empty sequence for $id\n" unless length $seq;
    die "Invalid non-ACGT reference tag for $id: $seq\n" unless $seq =~ /^[ACGT]+$/;
    if (!defined $$len_ref) {
        $$len_ref = length($seq);
    } elsif (length($seq) != $$len_ref) {
        die "Inconsistent tag lengths in $input: $id has " . length($seq) .
            " bp, expected $$len_ref bp\n";
    }
    print {$OUT} "$id $seq 100 100\n";
    $$count_ref++;
}

sub usage {
    my ($exit) = @_;
    print STDERR <<'USAGE';
Usage:
  perl mk_HQ_ref_codom.pl [-i ref_tag] [-o HQ_ref_codom]

Defaults:
  input  = ref_tag
  output = HQ_ref_codom

This script is enzyme-agnostic and automatically preserves the upstream tag
length. The trailing "100 100" fields are retained for format
compatibility; codom_calling.pl does not interpret them as parent information.
USAGE
    exit($exit // 0);
}
