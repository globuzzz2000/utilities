#!/usr/bin/env python3
"""
PLastZ-style alignment script
Splits sequences and runs LastZ alignments in parallel with proper MAF output
"""

import os
import subprocess
import shlex
import tempfile
import shutil
from multiprocessing import Pool
from Bio import SeqIO

def run_command(cmd):
    """Execute a shell command and return the exit code."""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Warning: Command failed with exit code {result.returncode}: {cmd}")
            if result.stderr:
                print(f"Error: {result.stderr}")
        return result.returncode
    except Exception as e:
        print(f"Error executing command: {cmd} - {str(e)}")
        return -1

def create_plastZ_jobs(ref_path, query_path, temp_dir, samtools_cmd="samtools", lastz_cmd="lastz", lastz_options=None):
    """
    Create extraction and alignment jobs following PLastZ approach.
    
    Returns:
        tuple: (extraction_jobs, alignment_jobs)
    """
    extract_jobs = []
    align_jobs = []
    extracted_sequences = []
    pairs_done = []
    
    print("Creating PLastZ-style alignment jobs...")
    
    # Get sequence IDs from both files
    ref_contigs = [seq.id for seq in SeqIO.parse(ref_path, "fasta")]
    query_contigs = [seq.id for seq in SeqIO.parse(query_path, "fasta")]
    
    print(f"Found {len(ref_contigs)} reference sequences and {len(query_contigs)} query sequences")
    
    # Create extraction and alignment jobs for each sequence pair
    for ref_id in ref_contigs:
        # Extract reference sequence if not already done
        if ref_id not in extracted_sequences:
            out_fa_path = os.path.join(temp_dir, f"{ref_id}.fa")
            ref_quoted = shlex.quote(ref_path)
            out_quoted = shlex.quote(out_fa_path)
            extract_jobs.append(f"{samtools_cmd} faidx {ref_quoted} {ref_id} > {out_quoted}")
            extracted_sequences.append(ref_id)
        
        for query_id in query_contigs:
            # Skip if we've already processed this pair
            if [ref_id, query_id] in pairs_done:
                continue
            
            # Mark pair as processed
            pairs_done.append([ref_id, query_id])
            
            # Extract query sequence if not already done
            if query_id not in extracted_sequences:
                out_fa_path = os.path.join(temp_dir, f"{query_id}.fa")
                query_quoted = shlex.quote(query_path)
                out_quoted = shlex.quote(out_fa_path)
                extract_jobs.append(f"{samtools_cmd} faidx {query_quoted} {query_id} > {out_quoted}")
                extracted_sequences.append(query_id)
            
            # Create alignment job
            tmp_output = os.path.join(temp_dir, f"{ref_id}_V_{query_id}.maf")  # Changed to .maf
            ref_fa = os.path.join(temp_dir, f"{ref_id}.fa")
            query_fa = os.path.join(temp_dir, f"{query_id}.fa")
            
            # Quote paths for shell safety
            ref_fa_quoted = shlex.quote(ref_fa)
            query_fa_quoted = shlex.quote(query_fa)
            tmp_output_quoted = shlex.quote(tmp_output)
            
            # Build LastZ command with PROPER MAF FORMAT
            base_cmd = f"{lastz_cmd} {ref_fa_quoted} {query_fa_quoted} --format=maf"
            
            # Add user options if provided
            if lastz_options:
                cmd = f"{base_cmd} {lastz_options} > {tmp_output_quoted}"
            else:
                cmd = f"{base_cmd} > {tmp_output_quoted}"
            
            align_jobs.append(cmd)
    
    print(f"Created {len(extract_jobs)} extraction jobs and {len(align_jobs)} alignment jobs")
    return extract_jobs, align_jobs

def run_plastZ_alignment(ref_path, query_path, output_path, processes=4, lastz_options=None, keep_temp=False, samtools_cmd="samtools", lastz_cmd="lastz"):
    """
    Run LastZ alignment using PLastZ approach with parallel processing.
    
    Args:
        ref_path (str): Path to reference FASTA file
        query_path (str): Path to query FASTA file
        output_path (str): Path for output MAF file
        processes (int): Number of parallel processes
        lastz_options (str): Additional LastZ options (--format=maf is automatically added)
        keep_temp (bool): Keep temporary files for debugging
    """
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix="plastZ_")
    print(f"Using temporary directory: {temp_dir}")
    
    try:
        # Create extraction and alignment jobs
        extract_jobs, align_jobs = create_plastZ_jobs(ref_path, query_path, temp_dir, samtools_cmd, lastz_cmd, lastz_options)
        
        # Run extraction jobs in parallel
        print(f"Extracting sequences using {processes} processes...")
        with Pool(processes=processes) as pool:
            extract_results = pool.map(run_command, extract_jobs)
        
        failed_extractions = sum(1 for r in extract_results if r != 0)
        if failed_extractions > 0:
            print(f"Warning: {failed_extractions}/{len(extract_jobs)} extraction jobs failed")
        
        # Run alignment jobs in parallel
        print(f"Running LastZ alignments using {processes} processes...")
        with Pool(processes=processes) as pool:
            align_results = pool.map(run_command, align_jobs)
        
        failed_alignments = sum(1 for r in align_results if r != 0)
        if failed_alignments > 0:
            print(f"Warning: {failed_alignments}/{len(align_jobs)} alignment jobs failed")
        
        # Combine all MAF files into final output with proper header
        print("Combining alignment results...")
        with open(output_path, 'w') as outfile:
            # Write MAF header
            outfile.write("##maf version=1 scoring=lastz.v1.04\n")
            outfile.write("# lastz alignment of sequences\n\n")
            
            maf_files = [f for f in os.listdir(temp_dir) if f.endswith('.maf')]
            for maf_file in sorted(maf_files):  # Sort for consistent output
                maf_path = os.path.join(temp_dir, maf_file)
                if os.path.exists(maf_path) and os.path.getsize(maf_path) > 0:
                    with open(maf_path, 'r') as infile:
                        content = infile.read()
                        # Skip header lines if they exist in individual files
                        lines = content.split('\n')
                        for line in lines:
                            if line.strip() and not line.startswith('##maf') and not line.startswith('# lastz'):
                                outfile.write(line + '\n')
        
        print(f"PLastZ alignment completed successfully!")
        print(f"Output saved to: {output_path}")
        
    except Exception as e:
        print(f"Error during PLastZ alignment: {str(e)}")
        raise
    
    finally:
        # Clean up temporary files
        if not keep_temp:
            print("Cleaning up temporary files...")
            shutil.rmtree(temp_dir)
            
            # Remove any .fai files created by samtools
            for fasta_file in [ref_path, query_path]:
                fai_file = f"{fasta_file}.fai"
                if os.path.exists(fai_file):
                    try:
                        os.remove(fai_file)
                    except:
                        pass
        else:
            print(f"Temporary files kept in: {temp_dir}")

if __name__ == "__main__":
    # Define your file paths
    reference_file = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Lyrata_filtered.fna"
    query_file = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Arenosa2_filtered.fna"
    output_file = "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Liftover/Lyrata_Arenosa2.maf"
    
    # Configuration
    num_processes = 6  # Adjust based on your CPU cores
    lastz_params = "--gapped --chain"  # Good defaults for genomic alignment
    keep_temp_files = False  # Set to True for debugging
    
    # Check if input files exist
    if not os.path.exists(reference_file):
        print(f"Error: Reference file not found: {reference_file}")
        exit(1)
        
    if not os.path.exists(query_file):
        print(f"Error: Query file not found: {query_file}")
        exit(1)
    
    # Check dependencies and find correct binary names
    def find_binary(names, test_args):
        """Try to find a binary by checking multiple possible names."""
        for name in names:
            try:
                result = subprocess.run([name] + test_args, capture_output=True)
                return name
            except FileNotFoundError:
                continue
        return None
    
    # Find samtools
    samtools_cmd = find_binary(["samtools"], ["--version"])
    if not samtools_cmd:
        print("Error: samtools not found.")
        exit(1)
    
    # Find lastz (lastz doesn't have --help, try with no arguments)
    lastz_cmd = find_binary(["lastz", "LASTZ", "LastZ"], [])
    if not lastz_cmd:
        print("Error: lastz not found.")
        exit(1)
    
    # Run the PLastZ alignment
    run_plastZ_alignment(reference_file, query_file, output_file, 
                        processes=num_processes, lastz_options=lastz_params, 
                        keep_temp=keep_temp_files, samtools_cmd=samtools_cmd, lastz_cmd=lastz_cmd)