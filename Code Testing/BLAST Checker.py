#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone BLAST script based on blast_processor.py
Blasts a list of sequences and outputs the two best hits for each.
"""

import os
import tempfile
import subprocess
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

class StandaloneBlastProcessor:
    """
    Standalone BLAST processor for short sequences.
    Based on the BlastProcessor class from blast_processor.py
    """
    
    # BLAST configuration (based on Config class from original)
    BLAST_WORD_SIZE = 7
    BLAST_EVALUE = 1000
    BLAST_REWARD = 2
    BLAST_PENALTY = -3
    BLAST_GAPOPEN = 5
    BLAST_GAPEXTEND = 2
    BLAST_MAX_TARGET_SEQS = 1000
    
    def __init__(self, db_path):
        """
        Initialize with BLAST database path.
        
        Args:
            db_path: Path to BLAST database
        """
        self.db_path = db_path
        self._temp_dir = None
        
    def _get_temp_dir(self):
        """Get or create temporary directory for BLAST files."""
        if self._temp_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="blast_temp_")
        return self._temp_dir
    
    def blast_short_seq(self, seq):
        """
        Run BLASTn for short sequences and return the two best e-values separately.
        
        Executes a BLASTn search optimized for short sequences and extracts
        the best and second-best e-values for specificity assessment. Only considers
        hits with 100% identity to ensure specificity.
        
        Args:
            seq: DNA sequence to BLAST against database
            
        Returns:
            Tuple of (best_evalue, second_best_evalue) or (None, None) if failed
        """
        if not seq or not isinstance(seq, str) or not seq.strip():
            logger.debug("Empty or invalid sequence provided to BLAST")
            return None, None

        tmp_filename = None
        try:
            # Create temporary file for query sequence
            temp_dir = self._get_temp_dir()
            
            with tempfile.NamedTemporaryFile(mode="w+", delete=False, dir=temp_dir) as tmp_query:
                tmp_query.write(f">seq\n{seq}\n")
                tmp_query.flush()
                tmp_filename = tmp_query.name

            # Execute BLASTn command with extended output format including identity
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
                    "-outfmt", "6 evalue pident qstart qend sstart send length"  # Extended format
                ],
                text=True,
                capture_output=True
            )
            
            if result.returncode != 0:
                error_msg = f"BLAST execution failed for sequence {seq[:20]}... (length: {len(seq)})"
                logger.error(error_msg)
                logger.error(f"BLAST stderr: {result.stderr}")
                return None, None

            # Parse BLAST output for e-values with 100% identity requirement
            try:
                valid_hits = []  # Store full hit information
                
                for line in result.stdout.strip().split("\n"):
                    if line.strip():
                        parts = line.strip().split("\t")
                        if len(parts) >= 7:  # evalue, pident, qstart, qend, sstart, send, length
                            evalue = float(parts[0])
                            pident = float(parts[1])
                            qstart = int(parts[2])
                            qend = int(parts[3])
                            sstart = int(parts[4])
                            send = int(parts[5])
                            length = int(parts[6])
                            
                            # Only include hits with 100% identity
                            if pident == 100.0:
                                hit_info = {
                                    'evalue': evalue,
                                    'pident': pident,
                                    'qstart': qstart,
                                    'qend': qend,
                                    'sstart': sstart,
                                    'send': send,
                                    'length': length
                                }
                                valid_hits.append(hit_info)
                            else:
                                logger.debug(f"Filtered BLAST hit with {pident}% identity (sequence {seq[:15]}...)")
                
                # Sort by e-value (best first)
                valid_hits = sorted(valid_hits, key=lambda x: x['evalue'])
                
            except ValueError as e:
                logger.warning(f"Error parsing BLAST output for sequence {seq[:20]}... (length: {len(seq)})")
                logger.debug(f"BLAST parsing error: {str(e)}")
                valid_hits = []

            if not valid_hits:
                logger.debug(f"No BLAST hits with 100% identity found for sequence {seq[:20]}... (length: {len(seq)})")
                return None, None, None, None

            # Extract best and second-best hits
            best_hit = valid_hits[0] if len(valid_hits) > 0 else None
            second_hit = valid_hits[1] if len(valid_hits) > 1 else None
            
            # Determine match type based on percent identity of best hit
            if best_hit:
                percent_identity = best_hit['pident']
                if percent_identity == 100.0:
                    match_type = 'blast_perfect'
                else:
                    match_type = 'blast_partial'  # This shouldn't happen since we filter for 100%
            else:
                match_type = None

            best_evalue = best_hit['evalue'] if best_hit else None
            second_evalue = second_hit['evalue'] if second_hit else None

            logger.debug(f"BLAST results for {seq[:20]}... -> Best: {best_evalue}, Second: {second_evalue}, Match type: {match_type}")
            
            return best_evalue, second_evalue, match_type, best_hit
            
        except Exception as e:
            error_msg = f"Unexpected BLAST error for sequence {seq[:20]}... (length: {len(seq)})"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}")
            return None, None, None, None
            
        finally:
            # Clean up temporary file
            if tmp_filename:
                try:
                    os.remove(tmp_filename)
                except OSError as e:
                    logger.debug(f"Failed to remove temp file {tmp_filename}: {e}")
    
    def cleanup(self):
        """Clean up temporary directory."""
        if self._temp_dir and os.path.exists(self._temp_dir):
            import shutil
            try:
                shutil.rmtree(self._temp_dir)
                logger.debug(f"Cleaned up temp directory: {self._temp_dir}")
            except OSError as e:
                logger.warning(f"Failed to clean up temp directory {self._temp_dir}: {e}")

def main():
    """Main function to run BLAST on the specified sequences."""
    
    # Sequences to BLAST
    sequences = [
        "CTAATGACCGCGGAAGAAG",
        "CGGGTGTGATTGGATCAAG", 
        "AGCGGCGCCACGGCTATCTTCCG",
        "CAACACCGAGTGATGAGAC",
        "CTGTTGTGGCTTCGGTTAG",
        "TGGCTCGACCGCTCTCGCTGCA",
        "GAACTGGGAAGTCAGATCAG",
        "GGATATTGGAGGACTTGACG",
        "CGCCCTTGGAGCTGCTCGTGGC",
        "CCGGAGTATCAAGCATGTC",
        "GTTTACCCTCTCGGCAATG",
        "CGTGGTTGCGCAACTCATCCCCGC",
        "GGAATCGGCGATCTACTTC",
        "CTTCAGAATCGGACTCCAC",
        "TGGAGAGGCGGCGGCTGAAGCT"
    ]
    
    # Database path from user_settings.json
    db_path = "/Users/jakob/.ddprimer/blast_db/Arabidopsis_thaliana"
    
    # Check if database exists
    if not os.path.exists(db_path + ".nin"):  # BLAST database files have extensions
        logger.error(f"BLAST database not found at: {db_path}")
        logger.error("Please check the database path and ensure BLAST database files exist.")
        return
    
    # Initialize BLAST processor
    blast_processor = StandaloneBlastProcessor(db_path)
    
    try:
        print("=" * 80)
        print("STANDALONE BLAST RESULTS")
        print("=" * 80)
        print(f"Database: {db_path}")
        print(f"Total sequences: {len(sequences)}")
        print("=" * 80)
        
        # Process each sequence
        for i, seq in enumerate(sequences, 1):
            print(f"\n[{i:2d}] Sequence: {seq}")
            print(f"     Length: {len(seq)} bp")
            
            # Run BLAST
            best_evalue, second_evalue, match_type, best_hit = blast_processor.blast_short_seq(seq)
            
            # Display results
            if best_evalue is not None and best_hit is not None:
                # Show match type and details
                print(f"     Match type:         {match_type}")
                print(f"     Percent identity:   {best_hit['pident']:.1f}%")
                print(f"     Best hit e-value:   {best_evalue:.2e}")
                print(f"     Alignment length:   {best_hit['length']} bp")
                print(f"     Query coverage:     {best_hit['qstart']}-{best_hit['qend']}")
                
                if second_evalue is not None:
                    print(f"     Second hit e-value: {second_evalue:.2e}")
                    ratio = best_evalue / second_evalue if second_evalue != 0 else float('inf')
                    print(f"     Specificity ratio:  {ratio:.2e}")
                    
                    # Interpret specificity
                    if ratio < 0.001:
                        print(f"     Specificity:        EXCELLENT (very specific)")
                    elif ratio < 0.01:
                        print(f"     Specificity:        GOOD (specific)")
                    elif ratio < 0.1:
                        print(f"     Specificity:        MODERATE (some cross-reactivity)")
                    else:
                        print(f"     Specificity:        POOR (high cross-reactivity)")
                else:
                    print(f"     Second hit e-value: None (only one hit)")
                    print(f"     Specificity ratio:  N/A (unique hit)")
                    print(f"     Specificity:        PERFECT (unique target)")
            else:
                print(f"     Match type:         None")
                print(f"     Best hit e-value:   None (no hits with 100% identity)")
                print(f"     Second hit e-value: None")
                print(f"     Specificity ratio:  N/A (no hits)")
                print(f"     Specificity:        NONE (no target found)")
        
        print("\n" + "=" * 80)
        print("BLAST ANALYSIS COMPLETE")
        print("=" * 80)
        
    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user.")
    except Exception as e:
        logger.error(f"Error during BLAST analysis: {str(e)}")
    finally:
        # Clean up
        blast_processor.cleanup()

if __name__ == "__main__":
    main()