"""
GC Content Analysis of Cancer-Related Genes

"""

import os

def read_fasta(filepath):
    """
    Read a FASTA file and return the DNA sequence as a string.
    Assumes a single-sequence FASTA file.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()

    sequence = ''
    for line in lines:
        if not line.startswith('>'):
            sequence += line.strip()
    return sequence.upper()


def compute_gc_content(seq):
    """
    Compute GC content percentage of a DNA sequence.
    """
    if not seq:
        return 0.0
    gc_count = seq.count('G') + seq.count('C')
    return round((gc_count / len(seq)) * 100, 2)


def analyze_gene_set(fasta_folder, gene_list):
    """
    Analyzes a list of genes and prints GC content metrics.
    """
    print("\n📊 GC Content Report for Selected Cancer Genes:\n")
    for gene in gene_list:
        file_path = os.path.join(fasta_folder, f"{gene}.fasta")
        if not os.path.exists(file_path):
            print(f"[!] Skipping {gene}: FASTA file not found.")
            continue

        sequence = read_fasta(file_path)
        gc_percent = compute_gc_content(sequence)
        total_len = len(sequence)
        at_count = sequence.count('A') + sequence.count('T')
        gc_count = sequence.count('G') + sequence.count('C')

        print(f"{gene}")
        print(f"    ➤ Length     : {total_len} bp")
        print(f"    ➤ GC Content : {gc_percent}%")
        print(f"    ➤ GC Count   : {gc_count}")
        print(f"    ➤ AT Count   : {at_count}\n")


if __name__ == "__main__":
    # Modify this list to include any other cancer genes you're analyzing
    genes_to_analyze = [
        "TP53",
        "BRCA1",
        "KRAS",
        "MYC",
        "FBXW7",
        "ARID1A",
        "PPP2R2A",
        "TET2"
    ]

    fasta_dir = "./fasta"
    analyze_gene_set(fasta_dir, genes_to_analyze)
