#!/usr/bin/env python3
"""Step 7: recover exact sample tags from serial-tag candidate-window FASTA files."""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from serial2brad16.common import find_pattern_hits, get_enzyme, iter_fasta, list_split_fasta, load_config, open_text, strip_fasta_suffix, write_fasta_record


def recover_one(inp: Path, out: Path, enzyme: dict, patterns, hit_policy: str):
    out.parent.mkdir(parents=True, exist_ok=True)
    total = no_hit = one_hit = multiple = emitted = 0
    with open_text(out, "wt") as handle:
        for header, seq in iter_fasta(inp):
            total += 1
            hits = [hit for hit in find_pattern_hits(seq, patterns) if len(hit[2]) == int(enzyme["tag_length"])]
            if not hits:
                no_hit += 1
                continue
            if len(hits) == 1:
                one_hit += 1
                selected = hits
            else:
                multiple += 1
                selected = [min(hits, key=lambda hit: (hit[3], hit[0], hit[1]))] if hit_policy == "first" else (hits if hit_policy == "all" else [])
            for hit_no, (start, end, tag, pattern_no) in enumerate(selected, 1):
                suffix = f"|tag0={start}-{end-1}|pattern={pattern_no}"
                if len(selected) > 1:
                    suffix += f"|hit={hit_no}"
                write_fasta_record(handle, header + suffix, tag)
                emitted += 1
    qc_path = out.parent / f"{strip_fasta_suffix(out.name)}.tag_recovery_qc.tsv"
    with open(qc_path, "w", encoding="utf-8") as handle:
        handle.write("input\tenzyme_id\tenzyme\ttag_length\tinput_windows\tno_hit\tone_hit\tmultiple_hits\temitted_tags\thit_policy\n")
        handle.write(f"{inp}\t{enzyme['id']}\t{enzyme['name']}\t{enzyme['tag_length']}\t{total}\t{no_hit}\t{one_hit}\t{multiple}\t{emitted}\t{hit_policy}\n")
    print(out)
    return inp, out, total, no_hit, one_hit, multiple, emitted


def main() -> int:
    parser = argparse.ArgumentParser(description="Serial2BRAD16 step 7: batch-recover exact enzyme-defined sample tags from all T1-T5 candidate FASTA files.")
    parser.add_argument("-d", "--input-dir", help="Directory mode (recommended: split5Data).")
    parser.add_argument("-r", "-i", "--input", dest="input", help="One T1-T5 candidate FASTA file.")
    parser.add_argument("-e", "--enzyme", required=True, type=int, choices=range(1, 17), help="One global enzyme preset ID (1-16 only).")
    parser.add_argument("-o", "--output", help="Output FASTA(.gz); valid only with one input file.")
    parser.add_argument("--outdir", default="tagData")
    parser.add_argument("--hit-policy", choices=("first", "unique", "all"), default="first")
    parser.add_argument("--config")
    args = parser.parse_args()

    if args.input_dir and args.input:
        parser.error("Use either directory mode (-d) or single-file mode (-r/-i), not both")
    if args.output and not args.input:
        parser.error("-o/--output is valid only with one input file")
    if args.input:
        inputs = [Path(args.input)]
        if not inputs[0].is_file():
            parser.error(f"Input does not exist: {inputs[0]}")
    else:
        try:
            inputs = list_split_fasta(Path(args.input_dir or "split5Data"))
        except ValueError as exc:
            parser.error(str(exc))
        if not inputs:
            parser.error("No _T1 through _T5 FASTA files were found")

    config = load_config(args.config)
    enzyme = get_enzyme(args.enzyme, config)
    patterns = [re.compile(pattern, re.IGNORECASE) for pattern in enzyme["patterns"]]
    outdir = Path(args.outdir)
    rows = []
    for inp in inputs:
        stem = strip_fasta_suffix(inp.name)
        out = Path(args.output) if args.output else outdir / f"{stem}.{enzyme['name']}.tags.fasta.gz"
        rows.append(recover_one(inp, out, enzyme, patterns, args.hit_policy))
    batch_dir = rows[0][1].parent if args.output else outdir
    batch_dir.mkdir(parents=True, exist_ok=True)
    with open(batch_dir / "tag_recovery_batch_qc.tsv", "w", encoding="utf-8") as handle:
        handle.write("input\toutput\tinput_windows\tno_hit\tone_hit\tmultiple_hits\temitted_tags\tenzyme_id\tenzyme\thit_policy\n")
        for inp, out, total, no_hit, one_hit, multiple, emitted in rows:
            handle.write(f"{inp}\t{out}\t{total}\t{no_hit}\t{one_hit}\t{multiple}\t{emitted}\t{enzyme['id']}\t{enzyme['name']}\t{args.hit_policy}\n")
    print(f"Recovered tags from {len(rows)} position file{'s' if len(rows) != 1 else ''}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
