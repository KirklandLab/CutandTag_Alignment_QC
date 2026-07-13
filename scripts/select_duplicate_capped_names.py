#!/usr/bin/env python3

import argparse
from collections import defaultdict

def is_flag_set(flag, bit):
    return (flag & bit) != 0

def main():
    parser = argparse.ArgumentParser(
        description="Select read names after capping duplicate paired-end fragments by exact fragment coordinates."
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--max-dups", type=int, required=True,
                        help="Maximum fragments to keep per identical fragment coordinate. Use <=0 to keep all.")
    parser.add_argument("--names-out", required=True)
    parser.add_argument("--metrics-out", required=True)
    args = parser.parse_args()

    max_dups = args.max_dups

    fragment_counts = defaultdict(int)
    kept_counts = defaultdict(int)

    total_fragments = 0
    kept_fragments = 0

    with open(args.names_out, "w") as names_out:
        import sys

        for line in sys.stdin:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            rnext = fields[6]
            pnext = int(fields[7])
            tlen = int(fields[8])

            # Use primary, mapped read1 records whose mate is also mapped.
            # This counts one fragment per paired-end template.
            if is_flag_set(flag, 0x4):      # read unmapped
                continue
            if is_flag_set(flag, 0x8):      # mate unmapped
                continue
            if is_flag_set(flag, 0x100):    # secondary
                continue
            if is_flag_set(flag, 0x800):    # supplementary
                continue
            if not is_flag_set(flag, 0x2):  # properly paired
                continue
            if not is_flag_set(flag, 0x40): # read1 only
                continue
            if tlen == 0:
                continue

            # Require same chromosome. RNEXT is "=" when mate maps to same ref.
            if not (rnext == "=" or rnext == rname):
                continue

            frag_start = min(pos, pnext)
            frag_end = frag_start + abs(tlen) - 1

            # Exact fragment coordinate duplicate key.
            # For Cut&Tag paired-end fragment analysis, this is the clearest practical definition.
            key = (rname, frag_start, frag_end)

            total_fragments += 1
            fragment_counts[key] += 1

            if max_dups <= 0 or kept_counts[key] < max_dups:
                kept_counts[key] += 1
                kept_fragments += 1
                names_out.write(qname + "\n")

    unique_fragment_positions = len(fragment_counts)
    duplicate_fragments_before = total_fragments - unique_fragment_positions
    removed_fragments = total_fragments - kept_fragments
    duplicate_fragments_after = kept_fragments - unique_fragment_positions

    pct_duplicate_before = (
        100 * duplicate_fragments_before / total_fragments
        if total_fragments > 0 else 0
    )

    pct_duplicate_after = (
        100 * duplicate_fragments_after / kept_fragments
        if kept_fragments > 0 else 0
    )

    with open(args.metrics_out, "w") as out:
        out.write(
            "sample\tmax_dups\ttotal_fragments_before\tunique_fragment_positions\t"
            "duplicate_fragments_before\tpct_duplicate_before\tkept_fragments_after\t"
            "removed_fragments\tduplicate_fragments_after\tpct_duplicate_after\n"
        )
        out.write(
            f"{args.sample}\t{max_dups}\t{total_fragments}\t{unique_fragment_positions}\t"
            f"{duplicate_fragments_before}\t{pct_duplicate_before:.4f}\t{kept_fragments}\t"
            f"{removed_fragments}\t{duplicate_fragments_after}\t{pct_duplicate_after:.4f}\n"
        )

if __name__ == "__main__":
    main()
