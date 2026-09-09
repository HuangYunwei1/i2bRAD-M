#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

# Extract Type IIB reference tags from a genome.
# The user selects one enzyme with -e (1-16). Recognition patterns and tag
# lengths are selected internally from -e.

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

my ($genome, $enzyme_selector, $outdir);
my $help = 0;
my $list_enzymes = 0;
GetOptions(
    'r=s'          => \$genome,
    'e=s'          => \$enzyme_selector,
    'o=s'          => \$outdir,
    'list-enzymes' => \$list_enzymes,
    'h|help'       => \$help,
) or usage(1);

if ($list_enzymes) {
    print_enzyme_table();
    exit 0;
}
usage(0) if $help;
usage(1, 'Required: -r <reference.fa[.gz]> -e <1-16> -o <output_directory>.')
    unless defined $genome && defined $enzyme_selector && defined $outdir;

my $enzyme_id = parse_enzyme_id($enzyme_selector);
my $spec = $ENZYMES{$enzyme_id};
my $enzyme = $spec->{name};
my $tag_len = $spec->{length};
my @patterns = map { qr/$_/ } @{ $spec->{patterns} };

make_path($outdir) unless -d $outdir;
my $ref_tag_path = "$outdir/ref_tag";
my $ref_path     = "$outdir/ref";
open(my $TAG, '>', $ref_tag_path) or die "Cannot write $ref_tag_path: $!\n";
open(my $REF, '>', $ref_path)     or die "Cannot write $ref_path: $!\n";
my $IN = open_fasta($genome);

my ($id, $seq) = (undef, '');
my $total = 0;
while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    if ($line =~ /^>(\S+)/) {
        process_record($id, $seq, \@patterns, $tag_len, $TAG, $REF, \$total) if defined $id;
        $id = $1;
        $seq = '';
    } else {
        die "Sequence before FASTA header in $genome\n" unless defined $id;
        $line =~ s/\s+//g;
        $seq .= uc($line);
    }
}
process_record($id, $seq, \@patterns, $tag_len, $TAG, $REF, \$total) if defined $id;

close $IN;
close $TAG;
close $REF;

print STDERR "Enzyme: $enzyme_id $enzyme (${tag_len} bp)\n";
print STDERR "Generated: $ref_tag_path\n";
print STDERR "Generated: $ref_path\n";
print STDERR "Reference tags: $total\n";

sub process_record {
    my ($id, $seq, $patterns, $length, $TAG, $REF, $total_ref) = @_;
    return unless defined $id && length($seq);
    my %hits;
    for my $pattern (@$patterns) {
        while ($seq =~ /($pattern)/g) {
            my $tag = $1;
            die "Internal enzyme definition error: extracted tag length " . length($tag) .
                " != expected $length\n" unless length($tag) == $length;
            my $end = pos($seq);  # 1-based end position
            $hits{$end} = $tag;
        }
    }

    my $count = 0;
    for my $end (sort { $a <=> $b } keys %hits) {
        $count++;
        $$total_ref++;
        my $header = ">$id-$count-$end";
        my $tag = $hits{$end};
        print {$TAG} "$header\n$tag\n";
        # Apply the SOAP2 padding design used by this workflow.
        print {$REF} "$header\n$tag", ('A' x 58), "\n";
    }
}

sub open_fasta {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        open($fh, '-|', 'gzip', '-dc', $path) or die "Cannot read $path via gzip: $!\n";
    } else {
        open($fh, '<', $path) or die "Cannot read $path: $!\n";
    }
    return $fh;
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
  perl Extract_cut_site.pl -r <reference.fa[.gz]> -e <1-16> -o <output_directory>

Required:
  -r   Reference genome FASTA (.gz supported)
  -e   One Type IIB enzyme ID (1-16)
  -o   Output directory (created if absent)

Optional:
  --list-enzymes  Show the 1-16 enzyme table
  -h              Show help

Outputs:
  <output_directory>/ref_tag   native reference tags
  <output_directory>/ref       padded reference used by SOAP2

The former -c and -l parameters have been removed. Recognition patterns and
correct tag lengths are selected internally from -e.
USAGE
    exit($exit // 0);
}
