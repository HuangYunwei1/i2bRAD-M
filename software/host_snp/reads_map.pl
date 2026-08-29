#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename);
use File::Path qw(make_path);

# Map sample-derived Type IIB tags to the padded reference with SOAP2.
# Historical -l is removed; tag length is selected internally from -e (1-16).

my %ENZYMES = (
     1 => { name => 'CspCI',  length => 33 },
     2 => { name => 'AloI',   length => 27 },
     3 => { name => 'BsaXI',  length => 27 },
     4 => { name => 'BaeI',   length => 28 },
     5 => { name => 'BcgI',   length => 32 },
     6 => { name => 'CjeI',   length => 28 },
     7 => { name => 'PpiI',   length => 27 },
     8 => { name => 'PsrI',   length => 27 },
     9 => { name => 'BplI',   length => 27 },
    10 => { name => 'FalI',   length => 27 },
    11 => { name => 'Bsp24I', length => 27 },
    12 => { name => 'HaeIV',  length => 27 },
    13 => { name => 'CjePI',  length => 27 },
    14 => { name => 'Hin4I',  length => 27 },
    15 => { name => 'AlfI',   length => 32 },
    16 => { name => 'BslFI',  length => 25 },
);

# Preserve the historical SOAP2 settings used by this workflow.
my $SOAP_MATCH_MODE = 4;
my $SOAP_MISMATCHES = 2;
my $SOAP_REPEAT_MODE = 0;

my $enzyme_selector;
my $help = 0;
my $list_enzymes = 0;
GetOptions(
    'e=s'          => \$enzyme_selector,
    'list-enzymes' => \$list_enzymes,
    'h|help'       => \$help,
) or usage(1);

if ($list_enzymes) {
    print_enzyme_table();
    exit 0;
}
usage(0) if $help;
usage(1, 'Required: -e <1-16>.') unless defined $enzyme_selector;

my $enzyme_id = parse_enzyme_id($enzyme_selector);
my $enzyme = $ENZYMES{$enzyme_id}{name};
my $tag_len = $ENZYMES{$enzyme_id}{length};

die "Required directory proc_data/ does not exist\n" unless -d 'proc_data';
die "Required SOAP2 reference soap/ref does not exist\n" unless -f 'soap/ref';
make_path('soap')          unless -d 'soap';
make_path('trash')         unless -d 'trash';
make_path('reads_mapping') unless -d 'reads_mapping';

check_executable('2bwt-builder');
check_executable('soap');
validate_reference('soap/ref', $tag_len);

run_cmd('2bwt-builder', 'soap/ref');

opendir(my $DIR, 'proc_data') or die "Cannot open proc_data/: $!\n";
my @files = sort grep {
    $_ ne '.' && $_ ne '..' &&
    -f "proc_data/$_" &&
    /\.(?:fa|fasta|fq|fastq)$/i
} readdir($DIR);
closedir $DIR;
die "No uncompressed FASTA/FASTQ sample files found in proc_data/\n" unless @files;

my @mapping_files;
for my $file (@files) {
    validate_sample_filename($file, $enzyme);
    my $sample = sample_name_from_file($file, $enzyme);
    my $soap_out = "soap/$sample.soap";
    my $unmapped = "trash/$sample.unmapped";

    run_cmd(
        'soap',
        '-a', "proc_data/$file",
        '-D', 'soap/ref.index',
        '-M', $SOAP_MATCH_MODE,
        '-r', $SOAP_REPEAT_MODE,
        '-u', $unmapped,
        '-v', $SOAP_MISMATCHES,
        '-o', $soap_out,
    );

    form_sample($sample, $soap_out, $tag_len);
    push @mapping_files, $sample;
}

# Keep the sample-column order explicit and tied to the files actually generated.
@mapping_files = sort @mapping_files;
open(my $ORDER, '>', 'sample_order.txt') or die "Cannot write sample_order.txt: $!\n";
print {$ORDER} "$_\n" for @mapping_files;
close $ORDER;

print STDERR "Enzyme: $enzyme_id $enzyme (${tag_len} bp)\n";
print STDERR "Processed ", scalar(@files), " sample file(s). Results: reads_mapping/\n";
print STDERR "Generated: sample_order.txt\n";

sub form_sample {
    my ($sample, $soap_file, $expected_len) = @_;

    my %id_to_tag;
    open(my $REF, '<', 'soap/ref') or die "Cannot read soap/ref: $!\n";
    my ($id, $seq) = (undef, '');
    while (my $line = <$REF>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        if ($line =~ /^>(\S+)/) {
            if (defined $id) {
                die "Reference sequence $id is shorter than expected tag length $expected_len\n"
                    if length($seq) < $expected_len;
                $id_to_tag{$id} = substr(uc($seq), 0, $expected_len);
            }
            $id = $1;
            $seq = '';
        } else {
            $line =~ s/\s+//g;
            $seq .= $line;
        }
    }
    if (defined $id) {
        die "Reference sequence $id is shorter than expected tag length $expected_len\n"
            if length($seq) < $expected_len;
        $id_to_tag{$id} = substr(uc($seq), 0, $expected_len);
    }
    close $REF;

    my (%reads, %depth);
    open(my $SOAP, '<', $soap_file) or die "Cannot read $soap_file: $!\n";
    while (my $line = <$SOAP>) {
        chomp $line;
        next unless length $line;
        my @f = split(/\t/, $line);
        die "Unexpected SOAP2 output format in $soap_file: fewer than 8 tab-separated columns\n" if @f < 8;
        my $ref_id = $f[7];
        next unless exists $id_to_tag{$ref_id};
        my $read_seq = uc($f[1] // '');
        die "SOAP2 query sequence length " . length($read_seq) .
            " does not match selected enzyme $enzyme expected tag length $expected_len in $soap_file\n"
            unless length($read_seq) == $expected_len;
        my $tag = $id_to_tag{$ref_id};
        $depth{$tag}++;
        push @{ $reads{$tag} }, $read_seq . ' ' . $f[0];
    }
    close $SOAP;

    my $outfile = "reads_mapping/$sample";
    open(my $OUT, '>', $outfile) or die "Cannot write $outfile: $!\n";
    for my $tag (sort keys %reads) {
        print {$OUT} "$_\n" for @{ $reads{$tag} };
        print {$OUT} "$tag   $depth{$tag}  ref\n";
    }
    close $OUT;
}

sub sample_name_from_file {
    my ($file, $selected_enzyme) = @_;
    my $name = basename($file);
    if ($name =~ /^(.*)\.\Q$selected_enzyme\E\.(?:fa|fasta|fq|fastq)$/i) {
        my $sample = $1;
        # Preserve the published/training BsaXI example: sample.BsaXI.fa -> sample.Bs.
        return "$sample.Bs" if lc($selected_enzyme) eq 'bsaxi';
        # For the other enzymes use an unambiguous enzyme-name suffix.
        return "$sample.$selected_enzyme";
    }
    $name =~ s/\.(?:fa|fasta|fq|fastq)$//i;
    return $name;
}

sub validate_sample_filename {
    my ($file, $selected_enzyme) = @_;
    my %known = map { lc($ENZYMES{$_}{name}) => $ENZYMES{$_}{name} } keys %ENZYMES;
    if ($file =~ /\.([A-Za-z0-9]+)\.(?:fa|fasta|fq|fastq)$/i) {
        my $label = lc($1);
        if (exists $known{$label} && $label ne lc($selected_enzyme)) {
            die "Sample file $file appears to contain $known{$label} tags, but reads_map.pl was run with $selected_enzyme\n";
        }
    }
}

sub validate_reference {
    my ($path, $expected_len) = @_;
    open(my $IN, '<', $path) or die "Cannot read $path: $!\n";
    my ($seen_header, $seen_seq) = (0, 0);
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        if ($line =~ /^>/) {
            $seen_header = 1;
        } else {
            die "$path contains sequence before FASTA header\n" unless $seen_header;
            my $seq = uc($line);
            $seq =~ s/\s+//g;
            die "$path contains a reference sequence shorter than selected tag length $expected_len\n"
                if length($seq) < $expected_len;
            my $tag = substr($seq, 0, $expected_len);
            die "$path reference tag prefix contains non-ACGT characters: $tag\n" unless $tag =~ /^[ACGT]+$/;
            $seen_seq++;
        }
    }
    close $IN;
    die "No reference sequences found in $path\n" unless $seen_seq;
}

sub check_executable {
    my ($exe) = @_;
    for my $dir (split(/:/, $ENV{PATH} // '')) {
        return if -x "$dir/$exe";
    }
    die "Cannot find required executable '$exe' in PATH\n";
}

sub run_cmd {
    my @cmd = @_;
    print STDERR 'RUN: ', join(' ', @cmd), "\n";
    system(@cmd);
    if ($? == -1) {
        die "Failed to execute $cmd[0]: $!\n";
    }
    if ($? & 127) {
        die sprintf("Command %s died with signal %d\n", $cmd[0], ($? & 127));
    }
    my $exit = $? >> 8;
    die "Command $cmd[0] failed with exit code $exit\n" if $exit != 0;
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
  From the work/ directory:
    perl reads_map.pl -e <1-16>

Expected workspace:
  proc_data/   uncompressed sample tag FASTA/FASTQ files
  soap/ref     padded reference generated by Extract_cut_site.pl with the same -e

Required:
  -e   One Type IIB enzyme ID (1-16). Tag length is selected internally.

Optional:
  --list-enzymes  Show the 1-16 enzyme table
  -h              Show help

Fixed SOAP2 settings retained from the workflow:
  -M 4   best-hit mode
  -v 2   maximum mismatches
  -r 0

Outputs:
  soap/            raw SOAP2 alignments / index files
  trash/           unmapped reads
  reads_mapping/   processed per-sample mapping files
  sample_order.txt exact mapping-file order used by downstream genotyping

BsaXI backward compatibility:
  sample.BsaXI.fa -> reads_mapping/sample.Bs
Other enzymes use an unambiguous suffix, e.g. sample.BcgI.fa -> sample.BcgI.
USAGE
    exit($exit // 0);
}
