#!/usr/bin/env python3

import os
from pathlib import Path

def parse_fasta(file_path):
    """Parse FASTA file and return list of (header, sequence) tuples with lengths."""
    sequences = []
    current_header = None
    current_sequence = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    seq_str = ''.join(current_sequence)
                    sequences.append((current_header, seq_str, len(seq_str)))
                
                # Start new sequence
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Don't forget the last sequence
        if current_header is not None:
            seq_str = ''.join(current_sequence)
            sequences.append((current_header, seq_str, len(seq_str)))
    
    return sequences

def filter_top_sequences(sequences, top_n=8):
    """Keep only the top N sequences by length."""
    # Sort by length (descending)
    sorted_sequences = sorted(sequences, key=lambda x: x[2], reverse=True)
    return sorted_sequences[:top_n]

def write_fasta(sequences, output_path):
    """Write sequences to FASTA file."""
    with open(output_path, 'w') as f:
        for header, sequence, length in sequences:
            f.write(f"{header}\n")
            # Write sequence in 80-character lines (standard FASTA format)
            for i in range(0, len(sequence), 80):
                f.write(f"{sequence[i:i+80]}\n")

def main():
    # Input file path
    input_path = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Lyrata.fna"
    
    # Create output file path (same directory, with "_top8" suffix)
    input_file = Path(input_path)
    output_path = input_file.parent / f"{input_file.stem}_top8{input_file.suffix}"
    
    print(f"Reading sequences from: {input_path}")
    
    # Check if input file exists
    if not os.path.exists(input_path):
        print(f"Error: Input file does not exist: {input_path}")
        return
    
    try:
        # Parse FASTA file
        sequences = parse_fasta(input_path)
        print(f"Found {len(sequences)} sequences")
        
        # Show sequence lengths
        print("\nAll sequence lengths:")
        for i, (header, seq, length) in enumerate(sequences, 1):
            print(f"{i:2d}. {header[:50]}{'...' if len(header) > 50 else ''} - {length:,} bp")
        
        # Filter to top 8
        top_sequences = filter_top_sequences(sequences, 8)
        
        print(f"\nKeeping top 8 sequences:")
        for i, (header, seq, length) in enumerate(top_sequences, 1):
            print(f"{i:2d}. {header[:50]}{'...' if len(header) > 50 else ''} - {length:,} bp")
        
        # Write filtered sequences
        write_fasta(top_sequences, output_path)
        print(f"\nFiltered FASTA written to: {output_path}")
        
        # Summary
        original_total = sum(seq[2] for seq in sequences)
        filtered_total = sum(seq[2] for seq in top_sequences)
        print(f"\nSummary:")
        print(f"Original: {len(sequences)} sequences, {original_total:,} total bp")
        print(f"Filtered: {len(top_sequences)} sequences, {filtered_total:,} total bp")
        print(f"Retained: {filtered_total/original_total*100:.1f}% of total sequence length")
        
    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    main()