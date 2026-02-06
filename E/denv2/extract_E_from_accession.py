#!/usr/bin/env python3
"""
Get COMPLETE E protein sequences from NCBI Protein database for DENV2.
Direct approach - avoids genome extraction issues.
"""

from Bio import Entrez, SeqIO
import time
import sys


def direct_protein_access(accession_ids, email):
    """
    Direct search for E proteins without genome extraction.
    """
    Entrez.email = email

    all_proteins = []

    # Search for DENV2 envelope proteins
    search_term = '(dengue virus 2[Organism] OR DENV-2[All Fields]) AND ' \
                  '(envelope[Title] OR glycoprotein E[Title] OR E protein[Title])'

    print("Searching for DENV2 envelope proteins in NCBI Protein database...")

    try:
        # First, search
        handle = Entrez.esearch(db="protein",term=search_term,retmax=10000,sort="relevance")
        results = Entrez.read(handle)
        handle.close()

        protein_count = len(results['IdList'])
        print(f"Found {protein_count} potential E protein entries")

        # Fetch in batches
        batch_size = 200
        for i in range(0, protein_count, batch_size):
            batch_ids = results['IdList'][i:i + batch_size]
            print(f"Fetching batch {i // batch_size + 1}/{(protein_count - 1) // batch_size + 1}")
            handle = Entrez.efetch(db="protein",id=",".join(batch_ids),rettype="fasta",retmode="text")
            batch_records = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            # Filter for complete E proteins
            for record in batch_records:
                # Check description
                desc = record.description.lower()

                # Must mention envelope/E protein
                if not ('envelope' in desc or ' e ' in desc or 'glycoprotein e' in desc):
                    continue

                # Must be reasonable length for E protein
                seq_len = len(record.seq)
                if not (450 <= seq_len <= 550):
                    continue

                # Check for common contaminants/issues
                seq_str = str(record.seq)
                if 'X' in seq_str:  # Skip sequences with unknown residues
                    continue

                if seq_str.count('M') < 2:  # Should have at least some methionines
                    continue

                all_proteins.append(record)

            time.sleep(0.5)

        print(f"\nRetrieved {len(all_proteins)} potential E proteins")

    except Exception as e:
        print(f"Error in direct search: {e}")

    return all_proteins


def deduplicate_proteins(protein_records):
    """Remove duplicate sequences based on sequence similarity."""
    unique_proteins = []
    seen_sequences = set()

    for record in protein_records:
        seq_str = str(record.seq)

        # Skip if we've seen this exact sequence
        if seq_str in seen_sequences:
            continue

        # Skip if sequence is too short
        if len(seq_str) < 450:
            continue

        seen_sequences.add(seq_str)
        unique_proteins.append(record)

    return unique_proteins


def main():
    # Configuration
    EMAIL = "maitysrithi@gmail.com"  # Already set to your email

    # Read accession IDs
    with open("DENV2_accession_ids.txt", "r") as f:
        accession_ids = [line.strip() for line in f if line.strip()]

    print(f"Loaded {len(accession_ids)} DENV2 accession IDs")

    # Direct protein search (only option)
    print("\n" + "=" * 60)
    print("DIRECT SEARCH FOR E PROTEINS IN PROTEIN DATABASE")
    print("=" * 60)

    e_proteins = direct_protein_access(accession_ids, EMAIL)

    if e_proteins:
        # Remove duplicates
        e_proteins = deduplicate_proteins(e_proteins)

        # Save all found proteins
        SeqIO.write(e_proteins, "DENV2_E_proteins_direct_search.fasta", "fasta")
        print(f"\nSaved {len(e_proteins)} E proteins to DENV2_E_proteins_direct_search.fasta")

        # Filter for complete ones (~495 aa)
        complete_e = [rec for rec in e_proteins if 480 <= len(rec.seq) <= 510]
        if complete_e:
            SeqIO.write(complete_e, "DENV2_E_proteins_complete.fasta", "fasta")
            print(f"Saved {len(complete_e)} complete E proteins (480-510 aa)")

        # Create a summary
        with open("DENV2_E_proteins_summary.txt", "w") as f:
            f.write(f"Total E proteins found: {len(e_proteins)}\n")
            f.write(f"Complete E proteins (480-510 aa): {len(complete_e)}\n\n")

            f.write("ALL E PROTEINS:\n")
            lengths = []
            for rec in e_proteins:
                lengths.append(len(rec.seq))
                f.write(f"{rec.id}: {len(rec.seq)} aa\n")

            if lengths:
                f.write(f"\nLENGTH STATISTICS:\n")
                f.write(f"Average: {sum(lengths) / len(lengths):.1f} aa\n")
                f.write(f"Min: {min(lengths)} aa\n")
                f.write(f"Max: {max(lengths)} aa\n")

        print(f"\nSTATISTICS:")
        print(f"  Total E proteins: {len(e_proteins)}")
        print(f"  Complete (480-510 aa): {len(complete_e)}")

        if lengths:
            print(f"  Average length: {sum(lengths) / len(lengths):.1f} aa")
            print(f"  Length range: {min(lengths)} - {max(lengths)} aa")

    else:
        print("\nNo E proteins found!")

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print("=" * 60)


if __name__ == "__main__":
    main()
