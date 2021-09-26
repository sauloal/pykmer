#!/usr/bin/env python3

import os
import gzip

ALFA = 'ACGT'

def iter_char(pos):
    for alfa in ALFA:
        if pos == 0:
            yield alfa
        else:
            for beta in iter_char(pos-1):
                yield alfa + beta

def create_sequence(kmer_len):
    for seq in iter_char(kmer_len-1):
        yield seq

def create_test(seq_name, kmer_len):
    fasta_file = f"{seq_name}-{kmer_len:02d}.fasta.gz"

    if os.path.exists(fasta_file): return

    with gzip.open(fasta_file, "wt") as fhd:
        for num, seq in enumerate(create_sequence(kmer_len)):
            fhd.write(f">{seq_name}-{kmer_len:02d}-{num+1:010d}\n{seq}\n")

def main():
    for kmer_len in [3, 5, 7, 9, 11, 13, 15, 17, 19, 21]:
        print(kmer_len)
        seq_name = "example"
        create_test(seq_name, kmer_len)

if __name__ == "__main__":
    main()
