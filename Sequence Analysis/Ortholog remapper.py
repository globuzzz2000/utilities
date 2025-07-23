#!/usr/bin/env python3
"""
Script to add Thaliana gene names from GFF file to ortholog mapping table
"""

import gzip
import csv
import re
from collections import defaultdict

def parse_gff_attributes(attributes_str):
    """Parse GFF attributes string to extract gene information"""
    
    gene_info = {}
    
    # Parse key=value; format (GFF3 style)
    for attr in attributes_str.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            gene_info[key.strip()] = value.strip().strip('"')
    
    return gene_info

def extract_gene_names_from_gff(gff_path):
    """Extract gene ID to gene name mapping from GFF file"""
    
    print(f"Reading GFF file: {gff_path}")
    gene_mapping = {}
    at_genes_found = 0
    
    try:
        with gzip.open(gff_path, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue
                
                # Split GFF line
                parts = line.split('\t')
                if len(parts) < 9:
                    continue
                
                feature = parts[2]
                attributes = parts[8]
                
                # Process gene features
                if feature.lower() == 'gene':
                    gene_info = parse_gff_attributes(attributes)
                    
                    # Extract AT gene ID and name based on the format we observed
                    gene_id = None
                    gene_name = None
                    
                    # Method 1: Check if ID contains AT gene (e.g., ID=gene-AT1G01010)
                    if 'ID' in gene_info and 'AT' in gene_info['ID'] and 'G' in gene_info['ID']:
                        id_field = gene_info['ID']
                        # Extract AT gene ID from ID=gene-AT1G01010 format
                        if id_field.startswith('gene-AT') and 'G' in id_field:
                            gene_id = id_field.replace('gene-', '')
                            at_genes_found += 1
                    
                    # Method 2: Check locus_tag field (e.g., locus_tag=AT1G01010)
                    elif 'locus_tag' in gene_info and gene_info['locus_tag'].startswith('AT') and 'G' in gene_info['locus_tag']:
                        gene_id = gene_info['locus_tag']
                        at_genes_found += 1
                    
                    # Get gene name from Name field
                    if gene_id and 'Name' in gene_info:
                        gene_name = gene_info['Name']
                        
                        # Store mapping
                        gene_mapping[gene_id] = gene_name
                        
                        # Debug: print first few mappings
                        if len(gene_mapping) <= 5:
                            print(f"  Found: {gene_id} -> {gene_name}")
                        
                # Show progress
                if line_num % 50000 == 0:
                    print(f"Processed {line_num} lines, found {len(gene_mapping)} AT gene mappings")
                    
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        return {}
    
    print(f"Found {len(gene_mapping)} AT gene name mappings")
    print(f"Total AT genes encountered: {at_genes_found}")
    return gene_mapping

def add_gene_names_to_mapping(csv_path, gene_mapping):
    """Add gene names to the ortholog mapping table"""
    
    print(f"Processing mapping table: {csv_path}")
    
    # Read original mapping
    rows = []
    try:
        with open(csv_path, 'r') as f:
            reader = csv.reader(f)
            rows = list(reader)
    except Exception as e:
        print(f"Error reading mapping table: {e}")
        return
    
    # Add gene names
    updated_rows = []
    found_count = 0
    
    for row in rows:
        if len(row) >= 2:
            species2_gene = row[0]
            thaliana_id = row[1]
            
            # Look up gene name
            gene_name = gene_mapping.get(thaliana_id, "")
            if gene_name:
                found_count += 1
            
            # Add gene name as third column
            updated_row = row + [gene_name]
            updated_rows.append(updated_row)
        else:
            updated_rows.append(row)
    
    # Write updated mapping
    output_path = csv_path.replace('.csv', '_with_gene_names.csv')
    
    try:
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(updated_rows)
        
        print(f"Updated mapping saved to: {output_path}")
        print(f"Found gene names for {found_count}/{len(rows)} entries")
        
    except Exception as e:
        print(f"Error writing output file: {e}")

def main():
    """Main function"""
    
    # File paths
    gff_path = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Thaliana.gff.gz"
    csv_path = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Ortholog_mapping.csv"
    
    # Extract gene names from GFF
    gene_mapping = extract_gene_names_from_gff(gff_path)
    
    if not gene_mapping:
        print("No gene mappings found. Please check the GFF file structure.")
        return
    
    # Show sample mappings
    print("\nSample gene mappings:")
    for i, (gene_id, gene_name) in enumerate(list(gene_mapping.items())[:5]):
        print(f"  {gene_id} -> {gene_name}")
    
    # Add gene names to mapping table
    add_gene_names_to_mapping(csv_path, gene_mapping)

if __name__ == "__main__":
    main()