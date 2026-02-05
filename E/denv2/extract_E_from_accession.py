#!/usr/bin/env python3
"""
Get COMPLETE E protein sequences from NCBI Protein database.
Direct approach - avoids genome extraction issues.
"""

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import time
import sys


def search_and_download_e_proteins(accession_ids, email, batch_size=20):
    """
    Search for E protein sequences in NCBI Protein database.
    """
    Entrez.email = email
    Entrez.api_key = None

    all_e_proteins = []
    failed_accessions = []

    print(f"Searching for E protein sequences for {len(accession_ids)} DENV2 genomes...")

    for i in range(0, len(accession_ids), batch_size):
        batch = accession_ids[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(accession_ids) - 1) // batch_size + 1

        print(f"\nBatch {batch_num}/{total_batches}: {len(batch)} genomes")

        try:
            # Search for protein sequences linked to these nucleotide accessions
            search_term = " OR ".join([f'{acc}[Nucleotide Accession]' for acc in batch])
            search_term += " AND dengue[Organism] AND envelope[Title]"

            # Search in protein database
            handle = Entrez.esearch(
                db="protein",
                term=search_term,
                retmax=batch_size * 5,  # Allow multiple proteins per genome
                retmode="xml"
            )

            search_results = Entrez.read(handle)
            handle.close()

            if not search_results['IdList']:
                print(f"  No E proteins found for this batch")
                continue

            # Fetch the protein sequences
            protein_ids = search_results['IdList']
            handle = Entrez.efetch(
                db="protein",
                id=",".join(protein_ids),
                rettype="fasta",
                retmode="text"
            )

            protein_records = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            # Filter for DENV2 E proteins
            batch_e_proteins = []
            for record in protein_records:
                desc = record.description.lower()

                # Check if it's DENV2 E protein
                is_denv2 = ('denv-2' in desc or 'dengue virus 2' in desc or
                            'dengue virus type 2' in desc or 'denv2' in desc)

                is_e_protein = ('envelope' in desc or 'glycoprotein e' in desc or
                                ' e protein' in desc)

                # Check length (E protein should be ~495 aa)
                length_ok = 450 <= len(record.seq) <= 550

                if is_denv2 and is_e_protein and length_ok:
                    batch_e_proteins.append(record)
                    print(f"  ✓ Found: {record.id} ({len(record.seq)} aa)")

            all_e_proteins.extend(batch_e_proteins)
            print(f"  Found {len(batch_e_proteins)} valid E proteins in this batch")

        except Exception as e:
            print(f"  Error: {e}")
            failed_accessions.extend(batch)

        # Respect NCBI rate limits
        if i + batch_size < len(accession_ids):
            time.sleep(1)

    return all_e_proteins, failed_accessions


def direct_protein_access(accession_ids, email):
    """
    Alternative: Direct search for E proteins without genome extraction.
    """
    Entrez.email = email

    all_proteins = []

    # Search for DENV2 envelope proteins
    search_term = '(dengue virus 2[Organism] OR DENV-2[All Fields]) AND ' \
                  '(envelope[Title] OR glycoprotein E[Title] OR E protein[Title])'

    print("Searching for DENV2 envelope proteins in NCBI Protein database...")

    try:
        # First, search
        handle = Entrez.esearch(
            db="protein",
            term=search_term,
            retmax=10000,  # Large enough
            sort="relevance"
        )

        results = Entrez.read(handle)
        handle.close()

        protein_count = len(results['IdList'])
        print(f"Found {protein_count} potential E protein entries")

        # Fetch in batches
        batch_size = 200
        for i in range(0, protein_count, batch_size):
            batch_ids = results['IdList'][i:i + batch_size]
            print(f"Fetching batch {i // batch_size + 1}/{(protein_count - 1) // batch_size + 1}")

            handle = Entrez.efetch(
                db="protein",
                id=",".join(batch_ids),
                rettype="fasta",
                retmode="text"
            )

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
    EMAIL = "maitysrithi@gmail.com"  # CHANGE THIS!

    # Read accession IDs
    with open("DENV2_accession_ids.txt", "r") as f:
        accession_ids = [line.strip() for line in f if line.strip()]

    print(f"Loaded {len(accession_ids)} DENV2 accession IDs")

    # OPTION 1: Direct protein search (recommended)
    print("\n" + "=" * 60)
    print("OPTION 1: Direct search for E proteins in Protein database")
    print("=" * 60)

    e_proteins = direct_protein_access(accession_ids, EMAIL)

    if e_proteins:
        # Remove duplicates
        e_proteins = deduplicate_proteins(e_proteins)

        # Save all found proteins
        SeqIO.write(e_proteins, "DENV2_E_proteins_direct_search.fasta", "fasta")
        print(f"\n✅ Saved {len(e_proteins)} E proteins to DENV2_E_proteins_direct_search.fasta")

        # Filter for complete ones (~495 aa)
        complete_e = [rec for rec in e_proteins if 480 <= len(rec.seq) <= 510]
        if complete_e:
            SeqIO.write(complete_e, "DENV2_E_proteins_complete.fasta", "fasta")
            print(f"✅ Saved {len(complete_e)} complete E proteins (480-510 aa)")

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

        print(f"\n📊 STATISTICS:")
        print(f"  Total E proteins: {len(e_proteins)}")
        print(f"  Complete (480-510 aa): {len(complete_e)}")

        if lengths:
            print(f"  Average length: {sum(lengths) / len(lengths):.1f} aa")
            print(f"  Length range: {min(lengths)} - {max(lengths)} aa")

    else:
        print("\n❌ No E proteins found!")

    # OPTION 2: Try linking to nucleotide accessions
    print("\n" + "=" * 60)
    print("OPTION 2: Search via nucleotide accession links")
    print("=" * 60)

    linked_e_proteins, failed = search_and_download_e_proteins(
        accession_ids[:50],  # Try first 50 as test
        EMAIL,
        batch_size=10
    )

    if linked_e_proteins:
        SeqIO.write(linked_e_proteins, "DENV2_E_proteins_linked.fasta", "fasta")
        print(f"\n✅ Saved {len(linked_e_proteins)} linked E proteins")


if __name__ == "__main__":
    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        print("Error: Biopython is required. Install with: pip install biopython")
        sys.exit(1)

    main()