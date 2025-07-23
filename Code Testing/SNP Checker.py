#!/usr/bin/env python3
"""
Combined BLAST and SNP checker script with strain name mapping.
First uses BLAST to find the best genomic location for each sequence,
then checks for SNPs at those locations and outputs genotype compatibility
with actual strain names instead of numeric IDs.
"""

from Bio import SeqIO
import pysam
import textwrap
import csv
import os
import tempfile
import subprocess
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# ---------- USER PATHS ----------
GENOME = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Thaliana.fna"
VCF    = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Thaliana.vcf.gz"
BLAST_DB = "/Users/jakob/.ddprimer/blast_db/Arabidopsis_thaliana"
OUTPUT_CSV = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/genotype_compatibility.csv"
ACCESSION_MAP = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Backup/Thaliana Accessions.csv"
# --------------------------------

# Sequences organized by genes
GENE_SEQUENCES = {
    "ILL6": [
        "CTAATGACCGCGGAAGAAG",
        "CGGGTGTGATTGGATCAAG", 
        "AGCGGCGCCACGGCTATCTTCCG"
    ],
    "ADC1": [
        "CAACACCGAGTGATGAGAC",
        "CTGTTGTGGCTTCGGTTAG",
        "TGGCTCGACCGCTCTCGCTGCA"
    ],
    "RLK902": [
        "GAACTGGGAAGTCAGATCAG",
        "GGATATTGGAGGACTTGACG",
        "CGCCCTTGGAGCTGCTCGTGGC"
    ],
    "FBH2": [
        "CCGGAGTATCAAGCATGTC",
        "GTTTACCCTCTCGGCAATG",
        "CGTGGTTGCGCAACTCATCCCCGC"
    ],
    "LON1": [
        "GGAATCGGCGATCTACTTC",
        "CTTCAGAATCGGACTCCAC",
        "TGGAGAGGCGGCGGCTGAAGCT"
    ]
}

class BlastProcessor:
    """BLAST processor for finding genomic locations."""
    
    # BLAST configuration
    BLAST_WORD_SIZE = 7
    BLAST_EVALUE = 1000
    BLAST_REWARD = 2
    BLAST_PENALTY = -3
    BLAST_GAPOPEN = 5
    BLAST_GAPEXTEND = 2
    BLAST_MAX_TARGET_SEQS = 1000
    
    def __init__(self, db_path):
        self.db_path = db_path
        self._temp_dir = None
        
    def _get_temp_dir(self):
        if self._temp_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="blast_temp_")
        return self._temp_dir
    
    def blast_sequence(self, seq):
        """
        Run BLASTn and return the best hit location.
        
        Returns:
            dict with 'blast_contig', 'blast_start', 'blast_end' or None if no hit
        """
        if not seq or not isinstance(seq, str) or not seq.strip():
            return None

        tmp_filename = None
        try:
            temp_dir = self._get_temp_dir()
            
            with tempfile.NamedTemporaryFile(mode="w+", delete=False, dir=temp_dir) as tmp_query:
                tmp_query.write(f">seq\n{seq}\n")
                tmp_query.flush()
                tmp_filename = tmp_query.name

            # Execute BLASTn with subject ID included in output
            result = subprocess.run(
                [
                    "blastn",
                    "-task", "blastn-short",
                    "-db", self.db_path,
                    "-query", tmp_filename,
                    "-word_size", str(self.BLAST_WORD_SIZE),
                    "-evalue", str(self.BLAST_EVALUE),
                    "-reward", str(self.BLAST_REWARD),
                    "-penalty", str(self.BLAST_PENALTY),
                    "-gapopen", str(self.BLAST_GAPOPEN),
                    "-gapextend", str(self.BLAST_GAPEXTEND),
                    "-max_target_seqs", str(self.BLAST_MAX_TARGET_SEQS),
                    "-outfmt", "6 sseqid evalue pident qstart qend sstart send length"
                ],
                text=True,
                capture_output=True
            )
            
            if result.returncode != 0:
                logger.error(f"BLAST execution failed for sequence {seq[:20]}...")
                return None

            # Parse BLAST output for best hit with 100% identity
            best_hit = None
            best_evalue = float('inf')
            
            for line in result.stdout.strip().split("\n"):
                if line.strip():
                    parts = line.strip().split("\t")
                    if len(parts) >= 8:
                        sseqid = parts[0]
                        evalue = float(parts[1])
                        pident = float(parts[2])
                        sstart = int(parts[5])
                        send = int(parts[6])
                        
                        # Only consider hits with 100% identity
                        if pident == 100.0 and evalue < best_evalue:
                            best_evalue = evalue
                            best_hit = {
                                'blast_contig': sseqid,
                                'blast_start': min(sstart, send) - 1,  # Convert to 0-based
                                'blast_end': max(sstart, send),        # End is exclusive
                                'evalue': evalue,
                                'strand': '+' if sstart < send else '-'
                            }
            
            return best_hit
            
        except Exception as e:
            logger.error(f"BLAST error for sequence {seq[:20]}...: {str(e)}")
            return None
            
        finally:
            if tmp_filename:
                try:
                    os.remove(tmp_filename)
                except OSError:
                    pass
    
    def cleanup(self):
        if self._temp_dir and os.path.exists(self._temp_dir):
            import shutil
            try:
                shutil.rmtree(self._temp_dir)
            except OSError:
                pass

def load_accession_mapping(accession_file):
    """
    Load accession mapping from CSV file.
    First column: numeric ID, third column: strain name
    
    Returns:
        dict: mapping from numeric ID to strain name
    """
    id_to_name = {}
    
    try:
        with open(accession_file, 'r') as f:
            csv_reader = csv.reader(f)
            for row in csv_reader:
                if len(row) >= 3:
                    numeric_id = row[0].strip()
                    strain_name = row[2].strip()
                    # Remove quotes if present
                    strain_name = strain_name.strip('"').strip("'")
                    id_to_name[numeric_id] = strain_name
        
        print(f"Loaded {len(id_to_name)} accession mappings")
        return id_to_name
        
    except FileNotFoundError:
        print(f"Warning: Accession mapping file not found: {accession_file}")
        print("Will use numeric IDs in output")
        return {}
    except Exception as e:
        print(f"Warning: Error loading accession mapping: {e}")
        print("Will use numeric IDs in output")
        return {}

def find_sequence_in_fasta(genome, query_seq, blast_coords):
    """
    Find the sequence in the local FASTA file using BLAST coordinates as reference.
    
    Args:
        genome: Dict of loaded FASTA sequences
        query_seq: The query sequence to find
        blast_coords: BLAST hit coordinates for reference
    
    Returns:
        Tuple of (contig, start, end) in local FASTA or None if not found
    """
    query_upper = query_seq.upper()
    
    # First try to find exact match in any contig
    for contig, record in genome.items():
        refseq = str(record.seq).upper()
        pos = refseq.find(query_upper)
        
        if pos != -1:
            start = pos
            end = pos + len(query_upper)
            return contig, start, end
    
    # If no exact match found, return None
    return None

def wrap(seq, width=100):
    return "\n".join(textwrap.wrap(seq, width))

def main():
    print("==== COMBINED BLAST + SNP CHECK RESULTS ====\n")

    # Check if BLAST database exists
    if not os.path.exists(BLAST_DB + ".nin"):
        logger.error(f"BLAST database not found at: {BLAST_DB}")
        return

    # Load accession mapping
    accession_map = load_accession_mapping(ACCESSION_MAP)
    
    # Load genome
    genome = SeqIO.to_dict(SeqIO.parse(GENOME, "fasta"))
    
    # Open VCF
    vcf = pysam.VariantFile(VCF)
    
    # Get all sample names from VCF
    sample_names = list(vcf.header.samples)
    print(f"Found {len(sample_names)} samples in VCF")
    
    # Initialize BLAST processor
    blast_processor = BlastProcessor(BLAST_DB)
    
    # Track which samples have variants for each gene
    gene_incompatible_samples = {gene: set() for gene in GENE_SEQUENCES.keys()}
    all_samples_with_variants = set()

    try:
        # Process each gene and its sequences
        for gene_name, sequences in GENE_SEQUENCES.items():
            print(f"\n=== Processing {gene_name} ===")
            
            for seq_idx, query in enumerate(sequences, start=1):
                query_upper = query.upper()
                
                print(f"\n{gene_name} Sequence {seq_idx}: {query}")
                print(f"Length: {len(query)} bp")
                
                # First, find best genomic location using BLAST
                blast_hit = blast_processor.blast_sequence(query)
                
                if not blast_hit:
                    print(f"  No BLAST hit found with 100% identity")
                    continue
                    
                blast_contig = blast_hit['blast_contig']
                blast_start = blast_hit['blast_start']
                blast_end = blast_hit['blast_end']
                evalue = blast_hit['evalue']
                strand = blast_hit['strand']
                
                print(f"  BLAST hit: {blast_contig}:{blast_start+1}-{blast_end} (strand: {strand}, e-value: {evalue:.2e})")

                # Now find the sequence in the local FASTA file
                fasta_location = find_sequence_in_fasta(genome, query, blast_hit)
                
                if not fasta_location:
                    print(f"  Sequence not found in local FASTA file")
                    continue
                    
                contig, start, end = fasta_location
                print(f"  Local FASTA location: {contig}:{start+1}-{end}")

                # Check if contig exists in VCF
                if contig not in vcf.header.contigs:
                    print(f"  Contig {contig} not present in VCF – skipping")
                    continue

                # Collect SNPs at the location
                snps = []
                for rec in vcf.fetch(contig, start, end):
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
                        
                        # Check each sample for variants at this position
                        for sample in rec.samples:
                            gt = rec.samples[sample]['GT']
                            if gt and any(allele != 0 for allele in gt if allele is not None):
                                gene_incompatible_samples[gene_name].add(sample)
                                all_samples_with_variants.add(sample)

                # Highlight SNPs in sequence
                seq_list = list(query_upper)
                for s in snps:
                    if 0 <= s["rel"] < len(seq_list):
                        seq_list[s["rel"]] = seq_list[s["rel"]].lower()
                highlighted = "".join(seq_list)

                # Report SNPs
                print(f"  Sequence with SNPs: {wrap(highlighted, 80)}")
                if snps:
                    print("  SNPs found:")
                    for s in snps:
                        print(f"    • pos_in_seq={s['rel']:<4}  ref={s['ref']}  alt={s['alt']}  AF={s['af']}  "
                              f"ID={s['id'] or '.'}  QUAL={s['qual']:.1f}  FILTER={s['filter']}")
                else:
                    print("  No SNPs found in this sequence.")

        # Summary by gene
        print(f"\n=== GENE COMPATIBILITY SUMMARY ===")
        for gene_name in GENE_SEQUENCES.keys():
            incompatible_count = len(gene_incompatible_samples[gene_name])
            compatible_count = len(sample_names) - incompatible_count
            print(f"{gene_name}: {compatible_count} compatible, {incompatible_count} incompatible")

        # Classify samples and convert to strain names
        compatible_samples = []
        incompatible_samples = []
        
        for sample in sample_names:
            strain_name = accession_map.get(sample, sample)  # Use strain name if available, otherwise use numeric ID
            
            if sample not in all_samples_with_variants:
                compatible_samples.append({
                    'strain': strain_name,
                    'sample_id': sample
                })
            else:
                incompatible_samples.append({
                    'strain': strain_name,
                    'sample_id': sample
                })

        print(f"\nOverall: {len(compatible_samples)} compatible, {len(incompatible_samples)} incompatible")

        # Write CSV with new format
        with open(OUTPUT_CSV, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Header row
            header = ['Accession'] + list(GENE_SEQUENCES.keys())
            writer.writerow(header)
            
            # Compatible samples first (no X marks)
            for sample_info in compatible_samples:
                row = [sample_info['strain']] + [''] * len(GENE_SEQUENCES)
                writer.writerow(row)
            
            # Incompatible samples with X marks for affected genes
            for sample_info in incompatible_samples:
                row = [sample_info['strain']]
                
                # Check each gene for incompatibility
                for gene_name in GENE_SEQUENCES.keys():
                    if sample_info['sample_id'] in gene_incompatible_samples[gene_name]:
                        row.append('X')
                    else:
                        row.append('')
                
                writer.writerow(row)

        print(f"\nResults written to: {OUTPUT_CSV}")
        
    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user.")
    except Exception as e:
        logger.error(f"Error during analysis: {str(e)}")
    finally:
        blast_processor.cleanup()

if __name__ == "__main__":
    main()