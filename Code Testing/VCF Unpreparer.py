#!/usr/bin/env python3
"""
VCF file modifier script to create test files that need normalization.

This script modifies a VCF file to:
1. Change chromosome names from accession numbers to standard format (chr1, chr2, etc.)
2. Introduce variants that need normalization (non-left-aligned indels, multiallelic variants)
3. Remove INFO/AF field to test AF field addition
4. Create issues that the file_preparator.py can detect and fix

Usage: python modify_vcf.py input.vcf output.vcf
"""

import sys
import re
import random
from pathlib import Path

def modify_vcf_for_testing(input_file, output_file):
    """
    Modify VCF file to create normalization and preparation issues.
    
    Args:
        input_file: Path to input VCF file
        output_file: Path to output modified VCF file
    """
    
    # Mapping from accession numbers to chromosome names
    chromosome_mapping = {
        'CP002684.1': 'chr1',
        'CP002685.1': 'chr2', 
        'CP002686.1': 'chr3',
        'CP002687.1': 'chr4',
        'CP002688.1': 'chr5'
    }
    
    # Counter for tracking modifications
    modifications = {
        'chromosomes_renamed': 0,
        'af_fields_removed': 0,
        'variants_made_non_normalized': 0,
        'multiallelic_created': 0
    }
    
    print(f"Reading from: {input_file}")
    print(f"Writing to: {output_file}")
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        variant_count = 0
        
        for line_num, line in enumerate(infile, 1):
            line = line.rstrip('\n')
            
            # Handle header lines
            if line.startswith('#'):
                if line.startswith('##contig='):
                    # Update contig lines with new chromosome names
                    for old_chr, new_chr in chromosome_mapping.items():
                        if old_chr in line:
                            # Extract length from original line
                            length_match = re.search(r'length=(\d+)', line)
                            if length_match:
                                length = length_match.group(1)
                                line = f"##contig=<ID={new_chr},length={length}>"
                            break
                elif line.startswith('##INFO=<ID=AF,'):
                    # Remove AF field definition to test AF field addition
                    modifications['af_fields_removed'] += 1
                    continue  # Skip this line
                
                outfile.write(line + '\n')
                continue
            
            # Handle variant lines
            fields = line.split('\t')
            if len(fields) < 8:
                outfile.write(line + '\n')
                continue
            
            variant_count += 1
            chrom, pos, var_id, ref, alt, qual, filter_field, info = fields[:8]
            
            # 1. Change chromosome names
            if chrom in chromosome_mapping:
                chrom = chromosome_mapping[chrom]
                modifications['chromosomes_renamed'] += 1
            
            # 2. Remove AF from INFO field to test AF field addition
            info_parts = info.split(';')
            info_parts = [part for part in info_parts if not part.startswith('AF=')]
            info = ';'.join(info_parts)
            
            # 3. Create normalization issues for some variants
            if variant_count % 7 == 0:  # Every 7th variant
                # Create non-left-aligned indel
                if len(ref) == 1 and len(alt) > 1:
                    # Make insertion non-left-aligned by adding context
                    ref = 'A' + ref  # Add leading base
                    alt = 'A' + alt
                    pos = str(int(pos) - 1)  # Adjust position
                    modifications['variants_made_non_normalized'] += 1
                elif len(ref) > 1 and len(alt) == 1:
                    # Make deletion non-left-aligned
                    ref = 'T' + ref
                    alt = 'T' + alt  
                    pos = str(int(pos) - 1)
                    modifications['variants_made_non_normalized'] += 1
            
            # 4. Create multiallelic variants occasionally
            if variant_count % 15 == 0 and len(ref) == 1 and len(alt) == 1:
                # Create multiallelic SNP
                bases = ['A', 'T', 'G', 'C']
                available_bases = [b for b in bases if b != ref and b != alt]
                if available_bases:
                    second_alt = random.choice(available_bases)
                    alt = f"{alt},{second_alt}"
                    modifications['multiallelic_created'] += 1
            
            # 5. Randomly introduce some problematic cases
            if variant_count % 20 == 0:
                # Add trailing whitespace (minor formatting issue)
                ref = ref + ' '
                alt = alt.strip()
            
            # Reconstruct the line
            new_fields = [chrom, pos, var_id, ref, alt, qual, filter_field, info]
            
            # Add remaining fields if they exist (FORMAT, samples, etc.)
            if len(fields) > 8:
                new_fields.extend(fields[8:])
            
            outfile.write('\t'.join(new_fields) + '\n')
    
    # Print summary of modifications
    print(f"\nModifications applied:")
    print(f"  - Chromosomes renamed: {modifications['chromosomes_renamed']}")
    print(f"  - AF field definitions removed: {modifications['af_fields_removed']}")
    print(f"  - Variants made non-normalized: {modifications['variants_made_non_normalized']}")
    print(f"  - Multiallelic variants created: {modifications['multiallelic_created']}")
    print(f"  - Total variants processed: {variant_count}")
    
    print(f"\nOutput file created: {output_file}")
    print("\nThis file now has issues that file_preparator.py can detect and fix:")
    print("  ✓ Non-standard chromosome names (CP002684.1 → chr1, etc.)")
    print("  ✓ Missing INFO/AF field")
    print("  ✓ Non-normalized variants (non-left-aligned indels)")
    print("  ✓ Multiallelic variants")
    print("  ✓ Will need bgzip compression and indexing")

def create_compatible_fasta_if_needed(original_vcf, output_fasta):
    """
    Optionally create a simple FASTA file with chromosome names that match the modified VCF.
    Only creates if user wants it.
    
    Args:
        original_vcf: Path to original VCF (to get chromosome info)
        output_fasta: Path to output FASTA file
    """
    
    response = input(f"\nDo you want to create a compatible FASTA file ({output_fasta})? [y/n]: ").strip().lower()
    if response not in ['y', 'yes']:
        print("Skipping FASTA file creation (using your existing FASTA file)")
        return
    
    # Read chromosome lengths from original VCF header
    chromosome_info = {}
    with open(original_vcf, 'r') as f:
        for line in f:
            if line.startswith('##contig='):
                # Extract ID and length
                id_match = re.search(r'ID=([^,]+)', line)
                length_match = re.search(r'length=(\d+)', line)
                if id_match and length_match:
                    chrom_id = id_match.group(1)
                    length = int(length_match.group(1))
                    chromosome_info[chrom_id] = length
            elif not line.startswith('#'):
                break
    
    # Mapping to new chromosome names
    chromosome_mapping = {
        'CP002684.1': 'chr1',
        'CP002685.1': 'chr2', 
        'CP002686.1': 'chr3',
        'CP002687.1': 'chr4',
        'CP002688.1': 'chr5'
    }
    
    print(f"Creating compatible FASTA file: {output_fasta}")
    
    with open(output_fasta, 'w') as f:
        for old_name, new_name in chromosome_mapping.items():
            if old_name in chromosome_info:
                length = chromosome_info[old_name]
                # Write FASTA header
                f.write(f">{new_name}\n")
                
                # Generate simple repetitive sequence (for testing purposes)
                # In real use, this would be the actual genomic sequence
                bases = 'ATGC'
                sequence = ''
                for i in range(length):
                    sequence += bases[i % 4]
                    if (i + 1) % 80 == 0:  # Line wrap at 80 characters
                        sequence += '\n'
                
                if not sequence.endswith('\n'):
                    sequence += '\n'
                
                f.write(sequence)
                print(f"  - {new_name}: {length:,} bp")
    
    print(f"FASTA file created with {len(chromosome_mapping)} chromosomes")

def main():
    """Main function to handle command line arguments and run the modification."""
    
    if len(sys.argv) != 3:
        print("Usage: python modify_vcf.py input.vcf output.vcf")
        print("\nThis script will:")
        print("  1. Change chromosome names to chr1, chr2, etc.")
        print("  2. Remove INFO/AF field")
        print("  3. Create non-normalized variants")
        print("  4. Create multiallelic variants")
        print("  5. Optionally generate a compatible FASTA file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Check if input file exists
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found!")
        sys.exit(1)
    
    # Create output directory if needed
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Modify the VCF file
    modify_vcf_for_testing(input_file, output_file)
    
    # Optionally create a compatible FASTA file
    fasta_file = output_path.with_suffix('.fasta')
    create_compatible_fasta_if_needed(input_file, str(fasta_file))
    
    print(f"\n" + "="*60)
    print("VCF FILE MODIFIED FOR TESTING:")
    print(f"  Modified VCF: {output_file}")
    print("\nYou can now test file_preparator.py with:")
    print(f"  - VCF file: {output_file}")
    print(f"  - Your existing FASTA file")
    print("\nThe preparator should detect and fix all the introduced issues.")

if __name__ == "__main__":
    main()