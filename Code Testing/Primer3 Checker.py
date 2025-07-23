#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone Primer3 Runner Script

This script runs Primer3 on a set of query sequences using the ddPrimer configuration.
It accesses the config from /Users/jakob/Applications/Git/ddPrimer/ddprimer/config/config.py
and processes the sequences through the Primer3Processor.
"""

import sys
import os
import logging
from typing import List, Dict

# Add the ddPrimer path to sys.path to import the config
sys.path.insert(0, '/Users/jakob/Applications/Git/ddPrimer')

try:
    # Import the ddPrimer config
    from ddprimer.config.config import Config
    # Import primer3 for direct usage
    import primer3
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Make sure you have primer3-py installed: pip install primer3-py")
    print("And that the ddPrimer path is correct")
    sys.exit(1)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Query sequences to test
QUERIES = [
    "CTAATGACCGCGGAAGAAGCTGCAAGAAGGAGATTAGCGGCGCCACGGCTATCTTCCGAAGTGATAACCGCCCGGAAAATTCCACATCCCGCGAGCAAAGGACCACTTCTTGATCCAATCACACCCG",
    "CAACACCGAGTGATGAGACACAATGGCTCGACCGCTCTCGCTGCAAATCACCGGATGCTTCACAGACTTCTGATCACAAACAAACCTAACCGAAGCCACAACAGC",
    "GAACTGGGAAGTCAGATCAGGCATCGCCCTTGGAGCTGCTCGTGGCTTAGACTATCTTCACTCACAAGACCCACTGAGCTCTCACGGAAACGTCAAGTCCTCCAATATCCT",
    "CCGGAGTATCAAGCATGTCGGATATGAACATGGAGAACCTTATGGAGGACTCTGTTGCTTTTAGGGTTCGGGCTAAACGTGGTTGCGCAACTCATCCCCGCAGCATTGCCGAGAGGGTAAAC",
    "GGAATCGGCGATCTACTTCGTTGGGTCATCGTGCTTTCTTCTGCTCTGAGCCTACCAATGGAGAGGCGGCGGCTGAAGCTGAAACTAAGGCGGTGGAGTCCGATTCTGAAGT"
]


class StandalonePrimer3Runner:
    """Standalone Primer3 runner using ddPrimer configuration."""
    
    def __init__(self):
        """Initialize with ddPrimer config."""
        self.config = Config.get_instance()
        logger.info("Initialized with ddPrimer configuration")
        
    def get_primer3_args(self) -> Dict:
        """Get Primer3 arguments from ddPrimer config."""
        return self.config.get_primer3_global_args()
    
    def run_primer3_on_sequence(self, sequence: str, sequence_id: str) -> Dict:
        """
        Run Primer3 on a single sequence.
        
        Args:
            sequence: DNA sequence to design primers for
            sequence_id: Identifier for the sequence
            
        Returns:
            Dictionary containing Primer3 results
        """
        logger.info(f"Designing primers for sequence {sequence_id} ({len(sequence)} bp)")
        
        # Prepare sequence arguments
        seq_args = {
            "SEQUENCE_ID": sequence_id,
            "SEQUENCE_TEMPLATE": sequence.upper(),
        }
        
        # Get global arguments from config
        global_args = self.get_primer3_args()
        
        try:
            # Run primer3
            result = primer3.bindings.design_primers(
                seq_args=seq_args,
                global_args=global_args
            )
            
            logger.info(f"Primer3 completed for {sequence_id}")
            return result
            
        except Exception as e:
            logger.error(f"Primer3 failed for {sequence_id}: {e}")
            return {}
    
    def format_primer_result(self, result: Dict, sequence_id: str) -> str:
        """
        Format Primer3 result for display.
        
        Args:
            result: Primer3 result dictionary
            sequence_id: Sequence identifier
            
        Returns:
            Formatted string representation of results
        """
        if not result:
            return f"\n=== {sequence_id} ===\nNo primers found or error occurred\n"
        
        output = [f"\n=== {sequence_id} ==="]
        
        # Check for errors
        if "PRIMER_ERROR" in result:
            output.append(f"PRIMER_ERROR: {result['PRIMER_ERROR']}")
        
        # Count primer pairs
        num_pairs = 0
        while f"PRIMER_LEFT_{num_pairs}_SEQUENCE" in result:
            num_pairs += 1
        
        output.append(f"Number of primer pairs returned: {num_pairs}")
        
        # Display each primer pair
        for i in range(num_pairs):
            output.append(f"\n--- Primer Pair {i+1} ---")
            
            # Pair penalty and product size
            if f"PRIMER_PAIR_{i}_PENALTY" in result:
                output.append(f"Pair Penalty: {result[f'PRIMER_PAIR_{i}_PENALTY']:.3f}")
            
            if f"PRIMER_PAIR_{i}_PRODUCT_SIZE" in result:
                output.append(f"Product Size: {result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']} bp")
            
            # Left primer
            if f"PRIMER_LEFT_{i}_SEQUENCE" in result:
                output.append(f"Forward Primer: {result[f'PRIMER_LEFT_{i}_SEQUENCE']}")
                
                if f"PRIMER_LEFT_{i}" in result:
                    start, length = result[f"PRIMER_LEFT_{i}"]
                    output.append(f"  Position: {start}-{start+length-1} (length: {length})")
                
                if f"PRIMER_LEFT_{i}_TM" in result:
                    output.append(f"  Tm: {result[f'PRIMER_LEFT_{i}_TM']:.1f}°C")
                
                if f"PRIMER_LEFT_{i}_PENALTY" in result:
                    output.append(f"  Penalty: {result[f'PRIMER_LEFT_{i}_PENALTY']:.3f}")
            
            # Right primer
            if f"PRIMER_RIGHT_{i}_SEQUENCE" in result:
                output.append(f"Reverse Primer: {result[f'PRIMER_RIGHT_{i}_SEQUENCE']}")
                
                if f"PRIMER_RIGHT_{i}" in result:
                    start, length = result[f"PRIMER_RIGHT_{i}"]
                    output.append(f"  Position: {start}-{start+length-1} (length: {length})")
                
                if f"PRIMER_RIGHT_{i}_TM" in result:
                    output.append(f"  Tm: {result[f'PRIMER_RIGHT_{i}_TM']:.1f}°C")
                
                if f"PRIMER_RIGHT_{i}_PENALTY" in result:
                    output.append(f"  Penalty: {result[f'PRIMER_RIGHT_{i}_PENALTY']:.3f}")
            
            # Internal oligo (probe) if present
            if f"PRIMER_INTERNAL_{i}_SEQUENCE" in result:
                output.append(f"Probe: {result[f'PRIMER_INTERNAL_{i}_SEQUENCE']}")
                
                if f"PRIMER_INTERNAL_{i}" in result:
                    start, length = result[f"PRIMER_INTERNAL_{i}"]
                    output.append(f"  Position: {start}-{start+length-1} (length: {length})")
                
                if f"PRIMER_INTERNAL_{i}_TM" in result:
                    output.append(f"  Tm: {result[f'PRIMER_INTERNAL_{i}_TM']:.1f}°C")
        
        return "\n".join(output)
    
    def display_config_summary(self):
        """Display a summary of the current configuration."""
        print("\n=== ddPrimer Configuration Summary ===")
        print(f"Primer Size Range: {self.config.PRIMER_MIN_SIZE}-{self.config.PRIMER_MAX_SIZE} bp")
        print(f"Primer Tm Range: {self.config.PRIMER_MIN_TM}-{self.config.PRIMER_MAX_TM}°C")
        print(f"Primer GC Range: {self.config.PRIMER_MIN_GC}-{self.config.PRIMER_MAX_GC}%")
        print(f"Product Size Range: {self.config.PRIMER_PRODUCT_SIZE_RANGE}")
        print(f"Max Primer Pairs: {self.config.MAX_PRIMER_PAIRS_PER_SEGMENT}")
        print("=" * 40)
    
    def run_all_queries(self) -> None:
        """Run Primer3 on all query sequences."""
        self.display_config_summary()
        
        logger.info(f"Processing {len(QUERIES)} query sequences...")
        
        for i, sequence in enumerate(QUERIES, 1):
            sequence_id = f"Query_{i}"
            
            # Run Primer3
            result = self.run_primer3_on_sequence(sequence, sequence_id)
            
            # Format and display results
            formatted_result = self.format_primer_result(result, sequence_id)
            print(formatted_result)
        
        logger.info("All sequences processed!")


def main():
    """Main function to run the standalone Primer3 runner."""
    print("Standalone Primer3 Runner using ddPrimer Configuration")
    print("=" * 55)
    
    try:
        # Create runner instance
        runner = StandalonePrimer3Runner()
        
        # Run on all queries
        runner.run_all_queries()
        
    except Exception as e:
        logger.error(f"Script failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()