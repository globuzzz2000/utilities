#!/usr/bin/env python3
"""
Locate exact-match query sequences in the Col-0 genome,
fetch overlapping SNPs from a bgzipped+indexed 1001-Genomes VCF,
and print the sequence with SNP sites in lowercase plus full SNP details.
"""

from Bio import SeqIO
import pysam
import textwrap

# ---------- USER PATHS ----------
GENOME = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Sequences/Col.fna"
VCF    = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Sequences/1001.vcf.gz"
# --------------------------------

QUERIES = [
    "CTAATGACCGCGGAAGAAGCTGCAAGAAGGAGATTAGCGGCGCCACGGCTATCTTCCGAAGTGATAACCGCCCGGAAAATTCCACATCCCGCGAGCAAAGGACCACTTCTTGATCCAATCACACCCG",
    "CAACACCGAGTGATGAGACACAATGGCTCGACCGCTCTCGCTGCAAATCACCGGATGCTTCACAGACTTCTGATCACAAACAAACCTAACCGAAGCCACAACAGC",
    "GAACTGGGAAGTCAGATCAGGCATCGCCCTTGGAGCTGCTCGTGGCTTAGACTATCTTCACTCACAAGACCCACTGAGCTCTCACGGAAACGTCAAGTCCTCCAATATCCT",
    "CCGGAGTATCAAGCATGTCGGATATGAACATGGAGAACCTTATGGAGGACTCTGTTGCTTTTAGGGTTCGGGCTAAACGTGGTTGCGCAACTCATCCCCGCAGCATTGCCGAGAGGGTAAAC",
    "GGAATCGGCGATCTACTTCGTTGGGTCATCGTGCTTTCTTCTGCTCTGAGCCTACCAATGGAGAGGCGGCGGCTGAAGCTGAAACTAAGGCGGTGGAGTCCGATTCTGAAGT"
]

# ----- mapping Col-0 accession IDs → VCF IDs -----
CONTIG_MAP = {
    "CP002684.1": "1",
    "CP002685.1": "2",
    "CP002686.1": "3",
    "CP002687.1": "4",
    "CP002688.1": "5",
    # chloroplast / mitochondrion are normally not in the nuclear-VCF:
    # "CP002689.1": "C",  # chloroplast
    # "CP002690.1": "M"   # mitochondrion
}
# -------------------------------------------------

def wrap(seq, width=100):
    return "\n".join(textwrap.wrap(seq, width))

def main():
    print("==== SNP CHECK RESULTS ====\n")

    # load genome once
    genome = SeqIO.to_dict(SeqIO.parse(GENOME, "fasta"))

    # open VCF (pysam auto-detects .tbi index)
    vcf = pysam.VariantFile(VCF)

    for idx, query in enumerate(QUERIES, start=1):
        query_upper = query.upper()
        hit_found = False

        for contig, record in genome.items():
            refseq = str(record.seq).upper()
            pos = refseq.find(query_upper)

            if pos == -1:
                continue  # not on this contig

            hit_found = True
            start, end = pos, pos + len(query_upper)

            # translate contig ID if necessary
            vcf_contig = CONTIG_MAP.get(contig, contig)
            if vcf_contig not in vcf.header.contigs:
                print(f"Sequence {idx}: contig {contig} ({vcf_contig}) not present in VCF – skipping.\n")
                break

            # collect SNPs
            snps = []
            for rec in vcf.fetch(vcf_contig, start, end):
                if rec.alts:
                    rel = rec.pos - start - 1  # 0-based within sequence
                    snps.append(
                        dict(
                            rel=rel,
                            ref=rec.ref,
                            alt=",".join(rec.alts),
                            id=rec.id,
                            af=rec.info.get("AF", "?"),
                            qual=rec.qual,
                            filter=";".join(rec.filter.keys()) if rec.filter else "PASS"
                        )
                    )

            # highlight SNPs
            seq_list = list(query_upper)
            for s in snps:
                if 0 <= s["rel"] < len(seq_list):
                    seq_list[s["rel"]] = seq_list[s["rel"]].lower()
            highlighted = "".join(seq_list)

            # report
            print(f"Sequence {idx}  (genome {contig}:{start+1}-{end},  VCF contig {vcf_contig})")
            print(wrap(highlighted))
            if snps:
                print("  SNPs:")
                for s in snps:
                    print(f"    • pos_in_seq={s['rel']:<4}  ref={s['ref']}  alt={s['alt']}  AF={s['af']}  "
                          f"ID={s['id'] or '.'}  QUAL={s['qual']:.1f}  FILTER={s['filter']}")
            else:
                print("  No SNPs found in this sequence.")
            print()
            break  # stop scanning other contigs for this query

        if not hit_found:
            print(f"Sequence {idx}: not found in the reference genome.\n")

if __name__ == "__main__":
    main()