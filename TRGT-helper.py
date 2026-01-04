#!/usr/bin/env python3

import argparse
import gzip
import sys
from pathlib import Path


def open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info_field(info_str):
    """Return dict of INFO key -> value (string)."""
    info = {}
    for field in info_str.split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        else:
            k, v = field, True
        info[k] = v
    return info


def main():
    parser = argparse.ArgumentParser(
        description="Extract TRGT allele sequences from a TRGT VCF into FASTA."
    )
    parser.add_argument("vcf", help="TRGT VCF (optionally gzipped)")
    parser.add_argument(
        "--prefix",
        default=None,
        help=(
            "Output prefix. If omitted: "
            "if 1 sample in VCF, use SAMPLE_TRID.fasta; "
            "else use TRID.fasta"
        ),
    )
    args = parser.parse_args()

    vcf_path = Path(args.vcf)

    samples = []
    records = []  # list of dicts: {sample, TRID, allele_index, seq}

    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header_fields = line.split("\t")
                samples = header_fields[9:]
                if len(samples) == 0:
                    sys.exit("No samples found in VCF FORMAT columns.")
                continue

            fields = line.split("\t")
            if len(fields) < 10:
                # not a proper genotype line
                continue

            chrom, pos, vid, ref, alts, qual, flt, info_str, fmt = fields[:9]
            sample_genotypes = fields[9:]

            info = parse_info_field(info_str)
            trid = info.get("TRID", None)
            if trid is None:
                # skip variants without TRID
                continue

            alt_list = alts.split(",")

            # FORMAT keys
            fmt_keys = fmt.split(":")

            for sname, geno in zip(samples, sample_genotypes):
                fmt_vals = geno.split(":")
                fmt_map = dict(zip(fmt_keys, fmt_vals))

                # TRGT puts motif-only region sequence in ALT alleles.
                # GT indexes refer to ALT alleles (1-based), 0 = REF.
                gt = fmt_map.get("GT", ".")
                if gt in {".", "./."}:
                    continue

                # split phased/unphased: 0/1, 1|2, etc.
                sep = "/" if "/" in gt else "|"
                allele_indices = gt.replace("|", "/").split("/")

                # For this use-case, just extract sequences directly from ALT,
                # ordered as allele-1, allele-2, etc.
                allele_seqs = []
                for idx_str in allele_indices:
                    try:
                        idx = int(idx_str)
                    except ValueError:
                        allele_seqs.append(None)
                        continue
                    if idx == 0:
                        # if 0 is used, you could choose ref or skip;
                        # here, use REF sequence as-is
                        allele_seqs.append(ref)
                    else:
                        # ALTs are 1-based
                        if 1 <= idx <= len(alt_list):
                            allele_seqs.append(alt_list[idx - 1])
                        else:
                            allele_seqs.append(None)

                # store a record per allele
                for i, seq in enumerate(allele_seqs, start=1):
                    if seq is None:
                        continue
                    records.append(
                        {
                            "sample": sname,
                            "TRID": trid,
                            "allele": i,
                            "seq": seq,
                        }
                    )

    if not records:
        sys.exit("No TRGT allele sequences found (no records with TRID).")

    # Determine output strategy
    unique_samples = sorted({r["sample"] for r in records})
    unique_trids = sorted({r["TRID"] for r in records})

    # Group records by TRID
    tr_to_records = {}
    for r in records:
        tr_to_records.setdefault(r["TRID"], []).append(r)

    for trid, recs in tr_to_records.items():
        if args.prefix is not None:
            # User-specified prefix: one FASTA per TRID using prefix + TRID
            out_prefix = f"{args.prefix}_{trid}"
        else:
            if len(unique_samples) == 1:
                # one sample; use SAMPLE_TRID.fasta
                out_prefix = f"{unique_samples[0]}_{trid}"
            else:
                # multiple samples; use TRID.fasta
                out_prefix = trid

        out_fasta = f"{out_prefix}.fasta"

        with open(out_fasta, "w") as out:
            for r in recs:
                # Header pattern:
                #   Sample_RFC1_Allele-1
                # or with TRID if multiple loci:
                #   Sample_TRID_Allele-1
                header = f"{r['sample']}_{r['TRID']}_Allele-{r['allele']}"
                out.write(f">{header}\n")

                # wrap sequence at 60 bp per line
                seq = r["seq"]
                for i in range(0, len(seq), 60):
                    out.write(seq[i : i + 60] + "\n")

        print(f"Wrote {out_fasta}", file=sys.stderr)


if __name__ == "__main__":
    main()

