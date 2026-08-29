#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use File::Copy qw(move);

# Multi-sample genotype calling for the Host-SNP workflow.
# Each sample is genotyped independently. Parent/progeny logic is absent.
# all_codom contains every attempted site. filter_of_all_codom retains sites
# satisfying BOTH the genotype call-rate threshold (-c) and MAF threshold (-m).

my $alpha      = 0.05;
my $min_depth  = 4;
my $max_depth  = 500;
my $error      = 0.01;
my $ref_file   = 'ref/HQ_ref_codom';
my $order_file = 'sample_order.txt';
my ($call_rate_threshold, $maf_threshold);
my $help = 0;

GetOptions(
    'c=f'          => \$call_rate_threshold,
    'm=f'          => \$maf_threshold,
    'a=f'          => \$alpha,
    'min-depth=i'  => \$min_depth,
    'max-depth=i'  => \$max_depth,
    'error=f'      => \$error,
    'ref=s'        => \$ref_file,
    'order=s'      => \$order_file,
    'h|help'       => \$help,
) or usage(1);
usage(0) if $help;
usage(1, 'Required: -c <genotype_call_rate> and -m <minor_allele_frequency>.')
    unless defined $call_rate_threshold && defined $maf_threshold;

die "-c must be a decimal between 0 and 1\n"
    unless $call_rate_threshold >= 0 && $call_rate_threshold <= 1;
die "-m must be a decimal between 0 and 0.5\n"
    unless $maf_threshold >= 0 && $maf_threshold <= 0.5;
die "Reference file $ref_file does not exist\n" unless -f $ref_file;
die "reads_mapping/ does not exist\n" unless -d 'reads_mapping';
die "--min-depth must be >= 1\n" unless $min_depth >= 1;
die "--max-depth must be greater than --min-depth\n" unless $max_depth > $min_depth;
die "--error must be > 0 and < 0.5\n" unless $error > 0 && $error < 0.5;

my $PP = ($alpha == 0.05) ? 3.84 : 5.90;  # preserve original two-threshold behavior

make_path('genotype') unless -d 'genotype';
my $all_file    = 'genotype/all_codom';
my $filter_file = 'genotype/filter_of_all_codom';
my $tmp_file    = 'genotype/.all_codom.tmp';

my @samples = load_sample_order($order_file);
die "No sample mapping files were found\n" unless @samples;

# Step 1: create Tag ID, reference tag, position, reference base backbone.
open(my $REF, '<', $ref_file) or die "Cannot read $ref_file: $!\n";
open(my $ALL, '>', $all_file) or die "Cannot write $all_file: $!\n";
my ($reference_tags, $reference_sites, $reference_tag_length) = (0, 0, undef);
while (my $line = <$REF>) {
    chomp $line;
    next unless $line =~ /\S/;
    my @f = split(/\s+/, $line);
    die "Malformed $ref_file line: $line\n" if @f < 2;
    my ($tag_id, $seq) = @f[0,1];
    $seq = uc($seq);
    die "Invalid reference tag sequence in $ref_file: $seq\n" unless $seq =~ /^[ACGT]+$/;
    if (!defined $reference_tag_length) {
        $reference_tag_length = length($seq);
    } elsif (length($seq) != $reference_tag_length) {
        die "Inconsistent reference tag lengths in $ref_file\n";
    }
    $reference_tags++;
    for my $pos (1 .. length($seq)) {
        my $base = substr($seq, $pos - 1, 1);
        print {$ALL} join(' ', $tag_id, $seq, $pos, $base), "\n";
        $reference_sites++;
    }
}
close $REF;
close $ALL;
die "No reference tags found in $ref_file\n" unless $reference_tags;

# Step 2: independently genotype every sample and append one column per sample.
for my $sample (@samples) {
    my $mapping = "reads_mapping/$sample";
    die "Sample listed in $order_file is missing: $mapping\n" unless -f $mapping;
    my $calls = genotype_sample($mapping, $min_depth, $max_depth, $error, $PP);

    open(my $IN,  '<', $all_file) or die "Cannot read $all_file: $!\n";
    open(my $TMP, '>', $tmp_file) or die "Cannot write $tmp_file: $!\n";
    while (my $line = <$IN>) {
        chomp $line;
        my @f = split(/\s+/, $line);
        my ($seq, $pos) = @f[1,2];
        my $gt = '--';
        if (exists $calls->{$seq} && exists $calls->{$seq}{$pos}) {
            $gt = $calls->{$seq}{$pos};
        }
        print {$TMP} "$line $gt\n";
    }
    close $IN;
    close $TMP;
    move($tmp_file, $all_file) or die "Cannot replace $all_file with temporary file: $!\n";
}

# Step 3: filter using BOTH call rate and MAF.
# call rate = number of samples with a valid diploid genotype / total samples.
# MAF = minimum NON-ZERO frequency among observed alleles in valid genotypes.
#       Monomorphic loci have MAF=0. Missing '--' samples do not contribute.
# Retention uses >= for both thresholds.
open(my $AIN, '<', $all_file) or die "Cannot read $all_file: $!\n";
open(my $FOUT, '>', $filter_file) or die "Cannot write $filter_file: $!\n";
my ($kept, $total) = (0, 0);
while (my $line = <$AIN>) {
    chomp $line;
    next unless $line =~ /\S/;
    $total++;
    my @f = split(/\s+/, $line);
    die "Malformed all_codom row: $line\n" if @f < 4 + @samples;
    my @gt = @f[4 .. $#f];
    die "Unexpected sample-column count in all_codom row\n" unless @gt == @samples;

    my ($call_rate, $maf, $called_samples) = calculate_call_rate_and_maf(\@gt);
    next if $called_samples == 0;  # MAF is not biologically defined without any genotype call.
    next if $call_rate < $call_rate_threshold;
    next if $maf < $maf_threshold;

    print {$FOUT} "$line\n";
    $kept++;
}
close $AIN;
close $FOUT;
unlink $tmp_file if -e $tmp_file;

print STDERR "Samples: ", scalar(@samples), "\n";
print STDERR "Reference tags: $reference_tags; tag length: $reference_tag_length bp; reference sites: $reference_sites\n";
print STDERR "Filter thresholds: call_rate >= $call_rate_threshold; MAF >= $maf_threshold\n";
print STDERR "Generated: $all_file\n";
print STDERR "Generated: $filter_file ($kept / $total sites retained)\n";

sub calculate_call_rate_and_maf {
    my ($gt_ref) = @_;
    my %allele_count = (A=>0, C=>0, G=>0, T=>0);
    my $called = 0;

    for my $gt (@$gt_ref) {
        next if !defined($gt) || $gt eq '--';
        die "Invalid genotype '$gt' in all_codom; expected -- or two A/C/G/T bases\n"
            unless $gt =~ /^[ACGT]{2}$/;
        $called++;
        my ($a1, $a2) = split //, $gt;
        $allele_count{$a1}++;
        $allele_count{$a2}++;
    }

    my $total_samples = scalar(@$gt_ref);
    my $call_rate = $total_samples ? $called / $total_samples : 0;
    return ($call_rate, 0, 0) if $called == 0;

    my $allele_total = 2 * $called;
    my @nonzero_freq = map { $allele_count{$_} / $allele_total }
                       grep { $allele_count{$_} > 0 } qw(A C G T);

    # A site with only one observed allele is monomorphic, so MAF=0.
    my $maf = 0;
    if (@nonzero_freq >= 2) {
        $maf = $nonzero_freq[0];
        for my $freq (@nonzero_freq[1 .. $#nonzero_freq]) {
            $maf = $freq if $freq < $maf;
        }
    }
    return ($call_rate, $maf, $called);
}

sub load_sample_order {
    my ($path) = @_;
    my @samples;

    if (-f $path) {
        open(my $IN, '<', $path) or die "Cannot read $path: $!\n";
        my %seen;
        while (my $line = <$IN>) {
            chomp $line;
            $line =~ s/^\s+|\s+$//g;
            next unless length $line;
            die "Duplicate sample mapping name in $path: $line\n" if $seen{$line}++;
            push @samples, $line;
        }
        close $IN;
        return @samples;
    }

    # Safe fallback if an older reads_map.pl did not create sample_order.txt.
    opendir(my $DIR, 'reads_mapping') or die "Cannot open reads_mapping/: $!\n";
    @samples = sort grep { $_ ne '.' && $_ ne '..' && -f "reads_mapping/$_" } readdir($DIR);
    closedir $DIR;
    open(my $OUT, '>', $path) or die "Cannot write $path: $!\n";
    print {$OUT} "$_\n" for @samples;
    close $OUT;
    print STDERR "Generated missing sample order file: $path\n";
    return @samples;
}

sub genotype_sample {
    my ($mapping, $min_d, $max_d, $e, $threshold) = @_;
    my %calls;
    my %counts;
    my $block_read_length;

    open(my $IN, '<', $mapping) or die "Cannot read $mapping: $!\n";
    while (my $line = <$IN>) {
        chomp $line;
        next unless $line =~ /\S/;
        my @f = split(/\s+/, $line);

        if (@f >= 3 && $f[-1] eq 'ref') {
            my $ref_seq = uc($f[0]);
            my $len = length($ref_seq);
            die "Invalid ref sequence in $mapping: $ref_seq\n" unless $ref_seq =~ /^[ACGT]+$/;
            if (defined $block_read_length && $block_read_length != $len) {
                die "Read/reference tag length mismatch in $mapping: reads=$block_read_length, ref=$len\n";
            }
            for my $pos (1 .. $len) {
                my $c = $counts{$pos} || { A=>0, C=>0, G=>0, T=>0 };
                $calls{$ref_seq}{$pos} = call_genotype($c, $min_d, $max_d, $e, $threshold);
            }
            %counts = ();
            undef $block_read_length;
        } else {
            my $read_seq = uc($f[0] // '');
            die "Invalid read sequence in $mapping: $read_seq\n" unless $read_seq =~ /^[ACGT]+$/;
            $block_read_length //= length($read_seq);
            die "Mixed read tag lengths within a mapping block in $mapping\n"
                if length($read_seq) != $block_read_length;
            for my $pos (1 .. length($read_seq)) {
                my $base = substr($read_seq, $pos - 1, 1);
                $counts{$pos}{$base}++;
            }
        }
    }
    close $IN;
    die "Mapping file $mapping ended before a terminating 'ref' line\n" if %counts;
    return \%calls;
}

sub call_genotype {
    my ($c, $min_d, $max_d, $e, $threshold) = @_;
    my @alleles = sort {
        (($c->{$b} // 0) <=> ($c->{$a} // 0)) || ($a cmp $b)
    } qw(A C G T);

    my ($a1, $a2) = @alleles[0,1];
    my ($n1, $n2) = (($c->{$a1} // 0), ($c->{$a2} // 0));
    my $depth = $n1 + $n2;
    return '--' if $depth < $min_d || $depth >= $max_d;

    my $hom_major = 1 - $e * 3 / 4;
    my $hom_minor = $e / 4;
    my $het_each  = 0.5 - $e / 4;
    my $lrt = 2 * (
        $n1 * log($hom_major) +
        $n2 * log($hom_minor) -
        $depth * log($het_each)
    );

    if ($lrt < 0 && -$lrt > $threshold) {
        return ($a1 le $a2) ? "$a1$a2" : "$a2$a1";
    }
    if ($lrt > $threshold) {
        return "$a1$a1";
    }
    return '--';
}

sub usage {
    my ($exit, $msg) = @_;
    print STDERR "$msg\n" if defined $msg;
    print STDERR <<'USAGE';
Usage:
  From the work/ directory:
    perl codom_calling.pl -c <call_rate> -m <MAF>

Required filter parameters (decimal):
  -c  Minimum genotype call rate, 0-1.
      call_rate = samples with genotype != '--' / total samples.
  -m  Minimum minor allele frequency, 0-0.5.
      MAF = the minimum NON-ZERO allele frequency among alleles actually observed
      in valid diploid genotype calls. Missing '--' samples are excluded.
      Monomorphic sites have MAF=0.

A site is written to filter_of_all_codom only when:
  call_rate >= -c  AND  MAF >= -m

Expected inputs:
  ref/HQ_ref_codom
  reads_mapping/<sample files>
  sample_order.txt              (normally generated by reads_map.pl)

Outputs:
  genotype/all_codom
  genotype/filter_of_all_codom

all_codom has no header:
  Col 1 = Tag ID
  Col 2 = Reference tag sequence
  Col 3 = Position within tag
  Col 4 = Reference base
  Col 5+ = Sample genotypes, in sample_order.txt order

Optional genotype-calling settings:
  -a              LRT significance setting [0.05 -> threshold 3.84; otherwise 5.90]
  --min-depth     Minimum major+second-allele depth for a call [4]
  --max-depth     Calls with depth >= this value are missing [500]
  --error         Sequencing error rate used by the LRT [0.01]
  --ref           Genotype reference [ref/HQ_ref_codom]
  --order         Sample order file [sample_order.txt]
  -h              Show help
USAGE
    exit($exit // 0);
}
