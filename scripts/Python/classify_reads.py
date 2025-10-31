#!/usr/bin/env python3
import sys, argparse, pysam
from itertools import groupby
from operator import attrgetter
from tqdm import tqdm

def nm_or_xm(rec):
    # Prefer NM; fall back to XM (Bowtie1)
    if rec.has_tag("NM"):
        return rec.get_tag("NM")
    elif rec.has_tag("XM"):
        return rec.get_tag("XM")
    else:
        # If neither: treat as unmapped-quality (push to worst)
        return 10**9

def groups_by_qname(bam):
    # pysam returns records sorted by QNAME after name-sort
    for qname, recs in tqdm(groupby(bam.fetch(until_eof=True), key=attrgetter("query_name"))):
        recs_list = [r for r in recs if not r.is_unmapped]
        yield qname, recs_list

def summarize_side(recs):
    """Return (min_mm, best_count, best_recs_subset) for a list of AlignedSegment."""
    if not recs:
        return (None, 0, [])
    mm_vals = [nm_or_xm(r) for r in recs]
    min_mm = min(mm_vals)
    best_recs = [r for r,mm in zip(recs, mm_vals) if mm == min_mm]
    return (min_mm, len(best_recs), best_recs)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mito", required=True)
    ap.add_argument("--decoy", required=True)
    ap.add_argument("--mito-out", required=True)
    ap.add_argument("--decoy-out", required=True)
    args = ap.parse_args()

    mito_bam   = pysam.AlignmentFile(args.mito, "rb")
    decoy_bam  = pysam.AlignmentFile(args.decoy, "rb")

    # Output templates: use mito header for mito_out, decoy header for decoy_out
    mito_out   = pysam.AlignmentFile(args.mito_out, "wb", template=mito_bam)
    decoy_out  = pysam.AlignmentFile(args.decoy_out, "wb", template=decoy_bam)

    # Iterate both name-sorted BAMs in lockstep by QNAME
    mito_iter  = groups_by_qname(mito_bam)
    decoy_iter = groups_by_qname(decoy_bam)

    try:
        m_name, m_recs = next(mito_iter)
    except StopIteration:
        m_name, m_recs = None, []
    try:
        d_name, d_recs = next(decoy_iter)
    except StopIteration:
        d_name, d_recs = None, []

    while m_name is not None or d_name is not None:
        # Align names
        if d_name is None or (m_name is not None and m_name < d_name):
            qname = m_name; m_cur = m_recs; d_cur = []
            try:
                m_name, m_recs = next(mito_iter)
            except StopIteration:
                m_name, m_recs = None, []
        elif m_name is None or d_name < m_name:
            qname = d_name; m_cur = []; d_cur = d_recs
            try:
                d_name, d_recs = next(decoy_iter)
            except StopIteration:
                d_name, d_recs = None, []
        else:
            qname = m_name  # equal
            m_cur = m_recs; d_cur = d_recs
            try:
                m_name, m_recs = next(mito_iter)
            except StopIteration:
                m_name, m_recs = None, []
            try:
                d_name, d_recs = next(decoy_iter)
            except StopIteration:
                d_name, d_recs = None, []

        # Summarize both sides
        m_min, m_best_n, m_best = summarize_side(m_cur)
        d_min, d_best_n, d_best = summarize_side(d_cur)

        # Both missing or unmapped -> drop
        if m_min is None and d_min is None:
            continue

        # Decision:
        # 1) Lower mismatches wins
        # 2) If equal mismatches: mito wins if mito_best_count < decoy_best_count
        # 3) Decoy is only kept if it wins by strictly fewer mismatches AND is unique at best (d_best_n == 1)
        if d_min is None or (m_min is not None and m_min < d_min):
            # mito wins
            for r in m_best:
                mito_out.write(r)
        elif m_min is None or d_min < m_min:
            # decoy strictly better; only keep if unique on decoy
            if d_best_n == 1:
                decoy_out.write(d_best[0])
        else:
            # m_min == d_min:
            if m_best_n < d_best_n:
                # mito wins tie by fewer best alignments
                for r in m_best:
                    mito_out.write(r)
            else:
                # decoy cannot win ties (by your rule); drop read
                pass

    mito_out.close(); decoy_out.close()
    mito_bam.close(); decoy_bam.close()

if __name__ == "__main__":
    main()
