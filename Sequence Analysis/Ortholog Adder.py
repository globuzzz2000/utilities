#!/usr/bin/env python3
"""
Script to update Arenosa GFF file with gene names from the mapping table
"""

import gzip
import csv
import re

def load_gene_name_mapping(csv_path):
    """Load the mapping from AARE701A_LOCUS ID to gene name"""
    
    print(f"Loading gene name mapping from: {csv_path}")
    mapping = {}
    
    try:
        with open(csv_path, 'r') as f:
            reader = csv.reader(f)
            for row_num, row in enumerate(reader, 1):
                if len(row) >= 3:
                    arenosa_id = row[0].strip()
                    gene_name = row[2].strip()
                    
                    if arenosa_id and gene_name:
                        mapping[arenosa_id] = gene_name
                    elif arenosa_id:
                        # Keep original ID if no gene name available
                        mapping[arenosa_id] = arenosa_id
                        
    except Exception as e:
        print(f"Error reading mapping file: {e}")
        return {}
    
    print(f"Loaded {len(mapping)} gene name mappings")
    return mapping

def parse_gff_attributes(attributes_str):
    """Parse GFF attributes string into a dictionary"""
    
    attributes = {}
    for attr in attributes_str.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key.strip()] = value.strip()
    
    return attributes

def rebuild_gff_attributes(attributes_dict):
    """Rebuild GFF attributes string from dictionary"""
    
    return ';'.join([f"{key}={value}" for key, value in attributes_dict.items()])

def update_arenosa_gff(input_gff_path, output_gff_path, gene_name_mapping):
    """Update the Arenosa GFF file with new gene names"""
    
    print(f"Updating GFF file: {input_gff_path}")
    print(f"Output file: {output_gff_path}")
    
    genes_updated = 0
    genes_not_found = 0
    total_lines = 0
    
    try:
        with gzip.open(input_gff_path, 'rt') as infile, \
             gzip.open(output_gff_path, 'wt') as outfile:
            
            for line_num, line in enumerate(infile, 1):
                total_lines += 1
                original_line = line
                line = line.strip()
                
                # Pass through comments and empty lines unchanged
                if line.startswith('#') or not line:
                    outfile.write(original_line)
                    continue
                
                # Split GFF line
                parts = line.split('\t')
                if len(parts) < 9:
                    outfile.write(original_line)
                    continue
                
                seqname, source, feature, start, end, score, strand, frame, attributes = parts
                
                # Process gene features with AARE701A_LOCUS
                if feature.lower() == 'gene' and 'AARE701A_LOCUS' in attributes:
                    
                    # Parse attributes
                    attr_dict = parse_gff_attributes(attributes)
                    
                    # Extract the AARE701A_LOCUS ID
                    locus_id = None
                    if 'locus_tag' in attr_dict:
                        locus_id = attr_dict['locus_tag']
                    elif 'ID' in attr_dict and 'AARE701A_LOCUS' in attr_dict['ID']:
                        # Extract from ID=gene-AARE701A_LOCUS1 format
                        id_field = attr_dict['ID']
                        if id_field.startswith('gene-'):
                            locus_id = id_field.replace('gene-', '')
                    
                    # Update the Name field if we have a mapping
                    if locus_id and locus_id in gene_name_mapping:
                        new_gene_name = gene_name_mapping[locus_id]
                        attr_dict['Name'] = new_gene_name
                        genes_updated += 1
                        
                        # Debug: print first few updates
                        if genes_updated <= 5:
                            print(f"  Updated: {locus_id} -> {new_gene_name}")
                            
                    elif locus_id:
                        genes_not_found += 1
                        if genes_not_found <= 5:
                            print(f"  No mapping found for: {locus_id}")
                    
                    # Rebuild the attributes string
                    new_attributes = rebuild_gff_attributes(attr_dict)
                    
                    # Write updated line
                    updated_parts = parts[:-1] + [new_attributes]
                    outfile.write('\t'.join(updated_parts) + '\n')
                    
                else:
                    # Write non-gene lines unchanged
                    outfile.write(original_line)
                
                # Progress reporting
                if line_num % 50000 == 0:
                    print(f"Processed {line_num} lines, updated {genes_updated} genes")
                    
    except Exception as e:
        print(f"Error processing GFF file: {e}")
        return
    
    print(f"\nUpdate complete!")
    print(f"Total lines processed: {total_lines}")
    print(f"Genes updated with new names: {genes_updated}")
    print(f"Genes not found in mapping: {genes_not_found}")
    print(f"Output written to: {output_gff_path}")

def main():
    """Main function"""
    
    # File paths
    mapping_csv_path = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Ortholog_mapping_with_gene_names.csv"
    input_gff_path = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Arenosa.gff.gz"
    output_gff_path = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Arenosa_with_gene_names.gff.gz"
    
    # Load gene name mapping
    gene_name_mapping = load_gene_name_mapping(mapping_csv_path)
    
    if not gene_name_mapping:
        print("No gene name mappings loaded. Exiting.")
        return
    
    # Show sample mappings
    print("\nSample gene name mappings:")
    for i, (locus_id, gene_name) in enumerate(list(gene_name_mapping.items())[:5]):
        print(f"  {locus_id} -> {gene_name}")
    
    # Update the GFF file
    update_arenosa_gff(input_gff_path, output_gff_path, gene_name_mapping)

if __name__ == "__main__":
    main()