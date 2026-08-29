#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

# Host-SNP / 2bRAD sample tag extraction.
# One Type IIB enzyme is selected per run with -e (1-16).
# Enzyme numbering/patterns follow the original 2bRADExtraction.pl and the
# MAP2B/2bRAD-M enzyme table used by this workflow.

my %ENZYMES = (
     1 => { name => 'CspCI',  length => 33, patterns => [ '[ACGT]{11}CAA[ACGT]{5}GTGG[ACGT]{10}', '[ACGT]{10}CCAC[ACGT]{5}TTG[ACGT]{11}' ] },
     2 => { name => 'AloI',   length => 27, patterns => [ '[ACGT]{7}GAAC[ACGT]{6}TCC[ACGT]{7}',   '[ACGT]{7}GGA[ACGT]{6}GTTC[ACGT]{7}' ] },
     3 => { name => 'BsaXI',  length => 27, patterns => [ '[ACGT]{9}AC[ACGT]{5}CTCC[ACGT]{7}',     '[ACGT]{7}GGAG[ACGT]{5}GT[ACGT]{9}' ] },
     4 => { name => 'BaeI',   length => 28, patterns => [ '[ACGT]{10}AC[ACGT]{4}GTA[CT]C[ACGT]{7}', '[ACGT]{7}G[AG]TAC[ACGT]{4}GT[ACGT]{10}' ] },
     5 => { name => 'BcgI',   length => 32, patterns => [ '[ACGT]{10}CGA[ACGT]{6}TGC[ACGT]{10}',   '[ACGT]{10}GCA[ACGT]{6}TCG[ACGT]{10}' ] },
     6 => { name => 'CjeI',   length => 28, patterns => [ '[ACGT]{8}CCA[ACGT]{6}GT[ACGT]{9}',      '[ACGT]{9}AC[ACGT]{6}TGG[ACGT]{8}' ] },
     7 => { name => 'PpiI',   length => 27, patterns => [ '[ACGT]{7}GAAC[ACGT]{5}CTC[ACGT]{8}',    '[ACGT]{8}GAG[ACGT]{5}GTTC[ACGT]{7}' ] },
     8 => { name => 'PsrI',   length => 27, patterns => [ '[ACGT]{7}GAAC[ACGT]{6}TAC[ACGT]{7}',    '[ACGT]{7}GTA[ACGT]{6}GTTC[ACGT]{7}' ] },
     9 => { name => 'BplI',   length => 27, patterns => [ '[ACGT]{8}GAG[ACGT]{5}CTC[ACGT]{8}' ] },
    10 => { name => 'FalI',   length => 27, patterns => [ '[ACGT]{8}AAG[ACGT]{5}CTT[ACGT]{8}' ] },
    11 => { name => 'Bsp24I', length => 27, patterns => [ '[ACGT]{8}GAC[ACGT]{6}TGG[ACGT]{7}',     '[ACGT]{7}CCA[ACGT]{6}GTC[ACGT]{8}' ] },
    12 => { name => 'HaeIV',  length => 27, patterns => [ '[ACGT]{7}GA[CT][ACGT]{5}[AG]TC[ACGT]{9}', '[ACGT]{9}GA[CT][ACGT]{5}[AG]TC[ACGT]{7}' ] },
    13 => { name => 'CjePI',  length => 27, patterns => [ '[ACGT]{7}CCA[ACGT]{7}TC[ACGT]{8}',      '[ACGT]{8}GA[ACGT]{7}TGG[ACGT]{7}' ] },
    14 => { name => 'Hin4I',  length => 27, patterns => [ '[ACGT]{8}GA[CT][ACGT]{5}[GAC]TC[ACGT]{8}', '[ACGT]{8}GA[CTG][ACGT]{5}[AG]TC[ACGT]{8}' ] },
    15 => { name => 'AlfI',   length => 32, patterns => [ '[ACGT]{10}GCA[ACGT]{6}TGC[ACGT]{10}' ] },
    16 => { name => 'BslFI',  length => 25, patterns => [ '[ACGT]{6}GGGAC[ACGT]{14}',               '[ACGT]{14}GTCCC[ACGT]{6}' ] },
);

my @input;
my ($enzyme_selector, $outdir, $outprefix);
my $q_control = 'yes';
my $format    = 'fa';
my $ncount    = 0.08;
my $quality   = 30;
my $percent   = 80;
my $qbase     = 33;
my $help      = 0;
my $list_enzymes = 0;

GetOptions(
    'i=s{1,2}'       => \@input,
    'e=s'            => \$enzyme_selector,
    'od=s'           => \$outdir,
    'op=s'           => \$outprefix,
    'qc=s'           => \$q_control,
    'fm=s'           => \$format,
    'n=f'            => \$ncount,
    'q=i'            => \$quality,
    'p=i'            => \$percent,
    'b=i'            => \$qbase,
    'list-enzymes'   => \$list_enzymes,
    'h|help'         => \$help,
) or usage(1);

if ($list_enzymes) {
    print_enzyme_table();
    exit 0;
}
usage(0) if $help;
usage(1, 'Required: -i <sample file>, -e <1-16>, -od <output dir>, -op <sample name>.')
    unless @input && defined $enzyme_selector && defined $outdir && defined $outprefix;

die "-qc must be yes or no\n" unless $q_control eq 'yes' || $q_control eq 'no';
die "-fm must be fa or fq\n" unless $format eq 'fa' || $format eq 'fq';
die "-p must be between 0 and 100\n" unless $percent >= 0 && $percent <= 100;
die "-n must be >= 0\n" unless $ncount >= 0;

my $enzyme_id = parse_enzyme_id($enzyme_selector);
my $spec = $ENZYMES{$enzyme_id};
my $enzyme = $spec->{name};
my $tag_len = $spec->{length};
my @patterns = map { qr/$_/ } @{ $spec->{patterns} };

make_path($outdir) unless -d $outdir;
my $outfile  = "$outdir/$outprefix.$enzyme.$format";
my $statfile = "$outdir/$outprefix.$enzyme.stat.xls";
open(my $OUT, '>', $outfile) or die "Cannot write $outfile: $!\n";

my ($input_reads, $enzyme_reads, $qc_reads) = (0, 0, 0);
my %input_types;

for my $file (@input) {
    my ($type, $IN) = detect_and_open($file);
    $input_types{$type}++;
    if ($type eq 'fastq') {
        process_fastq($IN, $file, $OUT, \@patterns, $tag_len,
                      \$input_reads, \$enzyme_reads, \$qc_reads);
    } else {
        die "-fm fq cannot be used with FASTA input ($file), because FASTA has no quality string\n"
            if $format eq 'fq';
        process_fasta($IN, $file, $OUT, \@patterns, $tag_len,
                      \$input_reads, \$enzyme_reads, \$qc_reads);
    }
    close $IN;
}
close $OUT;

open(my $STAT, '>', $statfile) or die "Cannot write $statfile: $!\n";
print {$STAT} "sample\tenzyme_id\tenzyme\ttag_length\tinput_reads_num\tenzyme_reads_num\tqc_reads_num\tpercent\n";
my $ratio = $input_reads ? sprintf('%.2f', $qc_reads / $input_reads * 100) : '0.00';
print {$STAT} join("\t", $outprefix, $enzyme_id, $enzyme, $tag_len,
                   $input_reads, $enzyme_reads, $qc_reads, "$ratio%"), "\n";
close $STAT;

print STDERR "Enzyme: $enzyme_id $enzyme (${tag_len} bp)\n";
print STDERR "Generated: $outfile\n";
print STDERR "Statistics: $statfile\n";

sub process_fastq {
    my ($IN, $file, $OUT, $patterns, $expected_len,
        $input_ref, $enzyme_ref, $qc_ref) = @_;
    while (1) {
        my $id = <$IN>;
        last unless defined $id;
        next if $id =~ /^\s*$/;
        my $seq  = <$IN>;
        my $plus = <$IN>;
        my $qual = <$IN>;
        die "Incomplete FASTQ record in $file\n" unless defined $seq && defined $plus && defined $qual;
        chomp($id, $seq, $plus, $qual);
        die "Malformed FASTQ header in $file: $id\n" unless $id =~ /^\@/;
        die "FASTQ sequence/quality length mismatch in $file for $id\n" unless length($seq) == length($qual);
        $$input_ref++;
        $seq = uc($seq);

        # Preserve the original type-3 behavior for long sequencing reads.
        if (length($seq) > 50) {
            $seq  = substr($seq,  0, 50);
            $qual = substr($qual, 0, 50);
        }

        my ($tag_seq, $tag_qual) = find_tag($seq, $qual, $patterns, $expected_len);
        next unless defined $tag_seq;
        $$enzyme_ref++;

        if ($q_control eq 'yes') {
            next unless pass_n_filter($tag_seq, $ncount);
            next unless pass_q_filter($tag_qual, $quality, $percent, $qbase);
        }
        $$qc_ref++;

        if ($format eq 'fa') {
            $id =~ s/^\@/>/;
            print {$OUT} "$id\n$tag_seq\n";
        } else {
            print {$OUT} "$id\n$tag_seq\n+\n$tag_qual\n";
        }
    }
}

sub process_fasta {
    my ($IN, $file, $OUT, $patterns, $expected_len,
        $input_ref, $enzyme_ref, $qc_ref) = @_;
    my ($id, $seq) = (undef, '');
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        if ($line =~ /^>(.*)$/) {
            process_fasta_record($id, $seq, $OUT, $patterns, $expected_len,
                                 $input_ref, $enzyme_ref, $qc_ref) if defined $id;
            $id = $1;
            $seq = '';
        } else {
            die "Sequence before FASTA header in $file\n" unless defined $id;
            $line =~ s/\s+//g;
            $seq .= uc($line);
        }
    }
    process_fasta_record($id, $seq, $OUT, $patterns, $expected_len,
                         $input_ref, $enzyme_ref, $qc_ref) if defined $id;
}

sub process_fasta_record {
    my ($id, $seq, $OUT, $patterns, $expected_len,
        $input_ref, $enzyme_ref, $qc_ref) = @_;
    $$input_ref++;
    $seq = substr($seq, 0, 50) if length($seq) > 50;
    my ($tag_seq) = find_tag($seq, undef, $patterns, $expected_len);
    return unless defined $tag_seq;
    $$enzyme_ref++;
    if ($q_control eq 'yes') {
        # FASTA has no quality values; only the N filter is meaningful.
        return unless pass_n_filter($tag_seq, $ncount);
    }
    $$qc_ref++;
    print {$OUT} ">$id\n$tag_seq\n";
}

sub find_tag {
    my ($seq, $qual, $patterns, $expected_len) = @_;
    for my $pattern (@$patterns) {
        if ($seq =~ /($pattern)/) {
            my $tag = $1;
            die "Internal enzyme definition error: extracted tag length " . length($tag) .
                " != expected $expected_len\n" unless length($tag) == $expected_len;
            my $start = $-[1];
            my $tag_qual = defined($qual) ? substr($qual, $start, $expected_len) : undef;
            return ($tag, $tag_qual);
        }
    }
    return;
}

sub detect_and_open {
    my ($file) = @_;
    my $fh = open_input($file);
    my $first;
    while (defined($first = <$fh>)) {
        next if $first =~ /^\s*$/;
        last;
    }
    die "Input file is empty: $file\n" unless defined $first;
    close $fh;
    my $type = ($first =~ /^>/) ? 'fasta' : ($first =~ /^\@/) ? 'fastq' : '';
    die "Cannot detect FASTA/FASTQ format for $file\n" unless $type;
    return ($type, open_input($file));
}

sub open_input {
    my ($file) = @_;
    my $fh;
    if ($file =~ /\.gz$/i) {
        open($fh, '-|', 'gzip', '-dc', $file) or die "Cannot read $file via gzip: $!\n";
    } else {
        open($fh, '<', $file) or die "Cannot read $file: $!\n";
    }
    return $fh;
}

sub pass_n_filter {
    my ($seq, $limit) = @_;
    my $n = ($seq =~ tr/N/N/);
    my $len = length($seq);
    return 0 if $len == 0;
    if ($limit > 0 && $limit < 1) {
        return ($n / $len <= $limit) ? 1 : 0;
    }
    return ($n <= $limit) ? 1 : 0;
}

sub pass_q_filter {
    my ($qual, $min_q, $min_percent, $base) = @_;
    return 0 unless defined $qual && length($qual);
    my $good = 0;
    for my $c (split //, $qual) {
        $good++ if ord($c) >= $min_q + $base;
    }
    return ($good >= length($qual) * $min_percent / 100) ? 1 : 0;
}

sub parse_enzyme_id {
    my ($raw) = @_;
    die "-e accepts exactly one enzyme ID (1-16) per Host-SNP run; multiple enzymes/17 are not supported here\n"
        unless defined $raw && $raw =~ /^\d+$/ && exists $ENZYMES{int($raw)};
    return int($raw);
}

sub print_enzyme_table {
    print "ID\tEnzyme\tTag_length\n";
    for my $id (sort { $a <=> $b } keys %ENZYMES) {
        print join("\t", $id, $ENZYMES{$id}{name}, $ENZYMES{$id}{length}), "\n";
    }
}

sub usage {
    my ($exit, $msg) = @_;
    print STDERR "$msg\n" if defined $msg;
    print STDERR <<'USAGE';
Usage:
  perl sample_extraction.pl -i <sample.fa/fq[.gz]> [-i accepts 1 or 2 files]
                            -e <1-16> -od <output_dir> -op <sample_name> [options]

Required:
  -i   Input FASTA/FASTQ path (.gz supported); one or two files accepted
  -e   One Type IIB enzyme ID (1-16); use --list-enzymes to display the table
  -od  Output directory (created if absent)
  -op  Output prefix / sample name

Optional:
  -qc  Quality control: yes|no [yes]
       FASTQ: N filter + quality filter; FASTA: N filter only
  -fm  Output format: fa|fq [fa]; fq requires FASTQ input
  -n   Maximum N ratio (0-1) or count [0.08]
  -q   Minimum base quality [30]
  -p   Minimum percent of bases reaching -q [80]
  -b   FASTQ quality ASCII base [33]
  --list-enzymes  Show the 1-16 enzyme table
  -h   Show help

Output:
  <output_dir>/<sample>.<Enzyme>.fa (or .fq)
  <output_dir>/<sample>.<Enzyme>.stat.xls

Note:
  One enzyme is analyzed per Host-SNP run. 17/AllEnzyme and comma-separated
  multi-enzyme selectors are intentionally not accepted in this pipeline.
USAGE
    exit($exit // 0);
}
