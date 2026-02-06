import os
from Bio import Entrez, SeqIO
import pandas as pd
from tqdm import tqdm
import time

# ==========================
# USER SETTINGS
# ==========================
Entrez.email = "maitysrithi@gmail.com"
Entrez.api_key = None

OUTPUT_DIR = "dengue_complete_sequences"

SEROTYPES = {
    "DENV1": "Dengue virus 1",
    "DENV2": "Dengue virus 2",
    "DENV3": "Dengue virus 3",
    "DENV4": "Dengue virus 4"
}

RETMAX = 5000  # batch size


# ==========================
# FUNCTIONS
# ==========================

def search_complete_genomes(virus_name):
    query = (
        f'"{virus_name}"[Organism] AND '
        'complete genome[Title] AND '
        'biomol_genomic[PROP]'
    )

    handle = Entrez.esearch(
        db="nucleotide",
        term=query,
        retmax=100000
    )
    record = Entrez.read(handle)
    handle.close()

    return record["IdList"]


def fetch_fasta_and_metadata(id_list, serotype):
    serotype_dir = os.path.join(OUTPUT_DIR, serotype)
    os.makedirs(serotype_dir, exist_ok=True)

    fasta_file = os.path.join(serotype_dir, f"{serotype}_complete.fasta")
    meta_file = os.path.join(serotype_dir, f"{serotype}_metadata.csv")
    accession_file = os.path.join(serotype_dir, f"{serotype}_accession_ids.txt")  # NEW

    metadata_rows = []
    accession_ids = []  # NEW: Store all accession IDs

    with open(fasta_file, "w") as fasta_out:
        for start in tqdm(range(0, len(id_list), RETMAX), desc=f"{serotype}"):
            batch_ids = id_list[start:start + RETMAX]
            ids = ",".join(batch_ids)

            # Fetch GenBank records
            handle = Entrez.efetch(
                db="nucleotide",
                id=ids,
                rettype="gb",
                retmode="text"
            )

            records = SeqIO.parse(handle, "genbank")

            for rec in records:
                # Write FASTA
                SeqIO.write(rec, fasta_out, "fasta")

                # Collect accession ID
                accession_ids.append(rec.id)  # NEW

                # Extract metadata
                meta = {
                    "accession": rec.id,
                    "length": len(rec.seq),
                    "description": rec.description,
                    "organism": rec.annotations.get("organism", ""),
                    "collection_date": "",
                    "country": "",
                    "host": ""
                }

                for feature in rec.features:
                    if feature.type == "source":
                        quals = feature.qualifiers
                        meta["collection_date"] = quals.get("collection_date", [""])[0]
                        meta["country"] = quals.get("country", [""])[0]
                        meta["host"] = quals.get("host", [""])[0]

                metadata_rows.append(meta)

            handle.close()
            time.sleep(0.4)

    # Save metadata
    df = pd.DataFrame(metadata_rows)
    df.to_csv(meta_file, index=False)

    # NEW: Save accession IDs to separate file
    with open(accession_file, "w") as acc_out:
        for acc_id in accession_ids:
            acc_out.write(f"{acc_id}\n")

    # Also create a combined accession file for all serotypes
    combined_accession_file = os.path.join(OUTPUT_DIR, "all_accession_ids.txt")
    with open(combined_accession_file, "a") as combined_out:
        for acc_id in accession_ids:
            combined_out.write(f"{serotype}\t{acc_id}\n")

    return len(accession_ids)


# ==========================
# MAIN LOOP (SEROTYPE ITERATION)
# ==========================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Create a summary file
    summary_file = os.path.join(OUTPUT_DIR, "download_summary.txt")

    with open(summary_file, "w") as summary:
        summary.write("Dengue Virus Complete Genome Download Summary\n")
        summary.write("=" * 50 + "\n")
        summary.write(f"Download Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    total_accessions = 0

    for serotype, virus_name in SEROTYPES.items():
        print(f"\n Processing {serotype}")

        id_list = search_complete_genomes(virus_name)
        print(f"  Found {len(id_list)} complete genomes")

        if len(id_list) == 0:
            print("  No complete genomes found")
            continue

        count = fetch_fasta_and_metadata(id_list, serotype)
        total_accessions += count

        # Update summary
        with open(summary_file, "a") as summary:
            summary.write(f"{serotype} ({virus_name}):\n")
            summary.write(f"  Total found: {len(id_list)}\n")
            summary.write(f"  Successfully downloaded: {count}\n")
            summary.write(f"  Accession list: {serotype}/{serotype}_accession_ids.txt\n")
            summary.write(f"  Metadata: {serotype}/{serotype}_metadata.csv\n")
            summary.write(f"  Sequences: {serotype}/{serotype}_complete.fasta\n")
            summary.write("-" * 30 + "\n")

        print(f"  Done: {serotype} ({count} sequences downloaded)")

    # Final summary
    with open(summary_file, "a") as summary:
        summary.write(f"\nTOTAL ACCESSIONS DOWNLOADED: {total_accessions}\n")
        summary.write(f"Combined accession list: all_accession_ids.txt\n")

    print(f"\n{'=' * 60}")
    print(f"DOWNLOAD COMPLETE")
    print(f"{'=' * 60}")
    print(f"Total sequences downloaded: {total_accessions}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Summary file: {os.path.join(OUTPUT_DIR, 'download_summary.txt')}")
    print(f"All accession IDs: {os.path.join(OUTPUT_DIR, 'all_accession_ids.txt')}")



if __name__ == "__main__":
    main()
