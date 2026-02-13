#!/usr/bin/env python3
"""
Verify and clean DENV4 E protein sequences.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import re


def verify_e_proteins(input_file, output_file):
    """
    Verify E protein sequences and save cleaned version.
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    print(f"Total sequences in input: {len(records)}")

    cleaned_records = []
    stats = {
        'total': 0,
        'complete': 0,
        'partial': 0,
        'unknown_residues': 0,
        'bad_length': 0
    }

    for record in records:
        stats['total'] += 1

        # Check for unknown residues
        seq_str = str(record.seq)
        if 'X' in seq_str or 'U' in seq_str or 'Z' in seq_str or 'B' in seq_str:
            stats['unknown_residues'] += 1
            continue

        # Check length
        seq_len = len(seq_str)
        if not (480 <= seq_len <= 510):
            stats['bad_length'] += 1
            continue

        # Check if description says "partial"
        desc_lower = record.description.lower()
        is_partial = 'partial' in desc_lower

        if is_partial:
            stats['partial'] += 1
            # Check if it might actually be complete despite label
            if seq_len >= 480:  # Close to full length
                # Remove "partial" from description
                new_desc = re.sub(r'\bpartial\b', '', record.description, flags=re.IGNORECASE)
                new_desc = re.sub(r'\s+', ' ', new_desc).strip()
                record.description = new_desc
                cleaned_records.append(record)
                stats['complete'] += 1
        else:
            stats['complete'] += 1
            cleaned_records.append(record)

    # Save cleaned sequences
    if cleaned_records:
        SeqIO.write(cleaned_records, output_file, "fasta")

    return stats, cleaned_records


def check_start_codons(records):
    """Check if sequences start with methionine."""
    with_met = 0
    without_met = 0

    for record in records:
        if str(record.seq).startswith('M'):
            with_met += 1
        else:
            without_met += 1

    return with_met, without_met


def main():
    input_file = "DENV4_E_proteins_complete.fasta"
    output_file = "DENV4_E_proteins_verified.fasta"

    print("Verifying DENV4 E protein sequences...")
    print("=" * 50)

    stats, cleaned_records = verify_e_proteins(input_file, output_file)

    print(f"\nSTATISTICS:")
    print(f"Total sequences processed: {stats['total']}")
    print(f"Complete sequences: {stats['complete']}")
    print(f"Marked as partial: {stats['partial']}")
    print(f"Skipped (unknown residues): {stats['unknown_residues']}")
    print(f"Skipped (bad length): {stats['bad_length']}")
    print(f"Saved to {output_file}: {len(cleaned_records)}")

    # Check start codons
    if cleaned_records:
        with_met, without_met = check_start_codons(cleaned_records)
        print(f"\nSTART CODON ANALYSIS:")
        print(f"Start with methionine (M): {with_met}")
        print(f"Don't start with M: {without_met}")

        # Show some examples
        print(f"\nSAMPLE SEQUENCES (first 3):")
        for i, rec in enumerate(cleaned_records[:3], 1):
            print(f"{i}. {rec.id}")
            print(f"   Length: {len(rec.seq)} aa")
            print(f"   Starts with: {rec.seq[0]}")
            print(f"   Description: {rec.description[:80]}...")


if __name__ == "__main__":
    main()