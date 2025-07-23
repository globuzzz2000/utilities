#!/usr/bin/env python3
"""
Simplified A. lyrata to A. arenosa Liftover Pipeline
For use when chain file and prepared VCF already exist
Optimized for Mac M1 with comprehensive quality validation
"""

# =============================================================================
# INPUT FILE CONFIGURATION - EDIT THESE PATHS
# =============================================================================

# Input files - Update these paths to match your actual files
INPUT_CONFIG = {
    # Your existing chain file (from previous MAF conversion)
    'chain_file': "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Liftover/lyrata_to_arenosa.chain",
    
    # Your prepared A. lyrata VCF file (can be .vcf or .vcf.gz)
    'vcf_file': "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Arenosa_Lyrata.vcf.gz",
    
    # Your A. arenosa reference genome (FASTA format)
    'target_fasta': "/Users/jakob/Library/Mobile Documents/com~apple~CloudDocs/Studium/Master/Master's Thesis/Primer Design/Arenosa.fna",
    
    # Output directory for all results
    'output_dir': "/Volumes/wimmer-nas-0327.synology.me/home/Drive/Liftover",
    
    # Prefix for output files
    'output_prefix': 'lyrata_to_arenosa',
    
    # Enable verbose logging (True/False)
    'verbose': True
}

# =============================================================================
# PIPELINE CODE - DO NOT MODIFY BELOW THIS LINE
# =============================================================================

import os
import sys
import gzip
import re
import subprocess
import logging
import argparse
from pathlib import Path
from datetime import datetime
from collections import defaultdict, Counter
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('liftover_pipeline.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class SimplifiedLiftoverPipeline:
    """Simplified liftover pipeline for pre-existing chain and VCF files"""
    
    def __init__(self, config):
        self.config = config
        self.validation_results = {}
        self.output_dir = Path(config['output_dir'])
        self.output_dir.mkdir(exist_ok=True)
        
        # File paths
        self.chain_file = Path(config['chain_file'])
        self.vcf_file = Path(config['vcf_file'])
        self.target_fasta = Path(config['target_fasta'])
        self.output_prefix = config['output_prefix']
        
        # Output files
        self.lifted_vcf = self.output_dir / f"{self.output_prefix}_lifted.vcf"
        self.unlifted_bed = self.output_dir / f"{self.output_prefix}_unlifted.bed"
        self.validation_dir = self.output_dir / "validation"
        self.validation_dir.mkdir(exist_ok=True)
        
    def validate_input_files(self):
        """Validate that all required input files exist and are properly formatted"""
        logger.info("=== Step 1: Input File Validation ===")
        
        validation = {
            'step': 'input_validation',
            'timestamp': datetime.now().isoformat(),
            'valid': True,
            'files': {},
            'warnings': [],
            'errors': []
        }
        
        # Validate chain file
        if not self.chain_file.exists():
            validation['errors'].append(f"Chain file not found: {self.chain_file}")
            validation['valid'] = False
        else:
            chain_stats = self._validate_chain_file_format()
            validation['files']['chain_file'] = {
                'path': str(self.chain_file),
                'exists': True,
                'size_mb': round(self.chain_file.stat().st_size / (1024*1024), 2),
                'format_valid': chain_stats['valid'],
                'chain_count': chain_stats.get('chain_count', 0)
            }
            if not chain_stats['valid']:
                validation['errors'].extend(chain_stats.get('errors', []))
                validation['valid'] = False
        
        # Validate VCF file
        if not self.vcf_file.exists():
            validation['errors'].append(f"VCF file not found: {self.vcf_file}")
            validation['valid'] = False
        else:
            vcf_stats = self._validate_vcf_file_format()
            validation['files']['vcf_file'] = {
                'path': str(self.vcf_file),
                'exists': True,
                'size_mb': round(self.vcf_file.stat().st_size / (1024*1024), 2),
                'format_valid': vcf_stats['valid'],
                'variant_count_sample': vcf_stats.get('variant_count', 0),
                'chromosomes': vcf_stats.get('chromosomes', [])
            }
            if not vcf_stats['valid']:
                validation['errors'].extend(vcf_stats.get('errors', []))
                validation['valid'] = False
        
        # Validate target FASTA
        if not self.target_fasta.exists():
            validation['errors'].append(f"Target FASTA not found: {self.target_fasta}")
            validation['valid'] = False
        else:
            fasta_stats = self._validate_fasta_file()
            validation['files']['target_fasta'] = {
                'path': str(self.target_fasta),
                'exists': True,
                'size_mb': round(self.target_fasta.stat().st_size / (1024*1024), 2),
                'format_valid': fasta_stats['valid'],
                'sequence_count': fasta_stats.get('sequence_count', 0),
                'sequences': fasta_stats.get('sequences', [])[:10]  # First 10 sequences
            }
            if not fasta_stats['valid']:
                validation['errors'].extend(fasta_stats.get('errors', []))
                validation['valid'] = False
        
        # Check for chromosome name compatibility
        if validation['valid']:
            compatibility_check = self._check_chromosome_compatibility(
                validation['files']['vcf_file'].get('chromosomes', []),
                validation['files']['target_fasta'].get('sequences', [])
            )
            validation['chromosome_compatibility'] = compatibility_check
            if compatibility_check['warnings']:
                validation['warnings'].extend(compatibility_check['warnings'])
        
        self._save_validation_report('input_validation', validation)
        
        if validation['valid']:
            logger.info("✓ All input files validated successfully")
            logger.info(f"  Chain file: {validation['files']['chain_file']['chain_count']} chains")
            logger.info(f"  VCF file: ~{validation['files']['vcf_file']['variant_count_sample']} variants (sample)")
            logger.info(f"  FASTA file: {validation['files']['target_fasta']['sequence_count']} sequences")
            if validation['warnings']:
                logger.warning(f"Warnings found: {len(validation['warnings'])}")
                for warning in validation['warnings']:
                    logger.warning(f"  - {warning}")
        else:
            logger.error("✗ Input validation failed")
            for error in validation['errors']:
                logger.error(f"  - {error}")
            raise ValueError("Input file validation failed")
        
        return validation
    
    def _validate_chain_file_format(self):
        """Validate chain file format"""
        stats = {'valid': False, 'chain_count': 0, 'errors': []}
        
        try:
            with open(self.chain_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith('chain '):
                        stats['chain_count'] += 1
                        parts = line.split()
                        
                        # Basic chain format validation
                        if len(parts) < 11:
                            stats['errors'].append(f"Line {line_num}: Invalid chain header format")
                            continue
                        
                        try:
                            # Validate numeric fields
                            score = int(parts[1])
                            t_size = int(parts[3])
                            t_start = int(parts[5])
                            t_end = int(parts[6])
                            q_size = int(parts[8])
                            q_start = int(parts[10])
                            q_end = int(parts[11]) if len(parts) > 11 else int(parts[10])
                        except ValueError:
                            stats['errors'].append(f"Line {line_num}: Invalid numeric values in chain header")
                    
                    # Sample first 1000 lines for validation
                    if line_num > 1000:
                        break
            
            stats['valid'] = stats['chain_count'] > 0 and len(stats['errors']) == 0
            
        except Exception as e:
            stats['errors'].append(f"Error reading chain file: {str(e)}")
        
        return stats
    
    def _validate_vcf_file_format(self):
        """Validate VCF file format"""
        stats = {'valid': False, 'variant_count': 0, 'chromosomes': set(), 'errors': []}
        
        try:
            with self._open_file(self.vcf_file) as f:
                has_header = False
                
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    if line.startswith('##'):
                        has_header = True
                        continue
                    elif line.startswith('#CHROM'):
                        continue
                    elif line:
                        stats['variant_count'] += 1
                        fields = line.split('\t')
                        
                        if len(fields) < 8:
                            stats['errors'].append(f"Line {line_num}: Insufficient fields in VCF line")
                            continue
                        
                        # Track chromosomes
                        chrom = fields[0]
                        stats['chromosomes'].add(chrom)
                    
                    # Sample first 1000 variants for validation
                    if stats['variant_count'] > 1000:
                        break
            
            stats['chromosomes'] = sorted(list(stats['chromosomes']))
            stats['valid'] = has_header and stats['variant_count'] > 0 and len(stats['errors']) == 0
            
        except Exception as e:
            stats['errors'].append(f"Error reading VCF file: {str(e)}")
        
        return stats
    
    def _validate_fasta_file(self):
        """Validate FASTA file format"""
        stats = {'valid': False, 'sequence_count': 0, 'sequences': [], 'errors': []}
        
        try:
            with self._open_file(self.target_fasta) as f:
                current_seq = None
                
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    if line.startswith('>'):
                        if current_seq:
                            stats['sequences'].append(current_seq)
                        current_seq = line[1:].split()[0]  # Get sequence ID
                        stats['sequence_count'] += 1
                    elif current_seq and line:
                        # Validate DNA sequence
                        if not re.match(r'^[ACGTNacgtn]+$', line):
                            stats['errors'].append(f"Line {line_num}: Invalid DNA sequence characters")
                    
                    # Sample first 100 sequences for validation
                    if stats['sequence_count'] > 100:
                        break
                
                if current_seq:
                    stats['sequences'].append(current_seq)
            
            stats['valid'] = stats['sequence_count'] > 0 and len(stats['errors']) == 0
            
        except Exception as e:
            stats['errors'].append(f"Error reading FASTA file: {str(e)}")
        
        return stats
    
    def _check_chromosome_compatibility(self, vcf_chroms, fasta_seqs):
        """Check compatibility between VCF chromosomes and FASTA sequences"""
        compatibility = {
            'vcf_chromosomes': vcf_chroms[:10],  # First 10
            'fasta_sequences': fasta_seqs[:10],  # First 10
            'potential_matches': 0,
            'warnings': []
        }
        
        # Simple matching check
        vcf_set = set(vcf_chroms)
        fasta_set = set(fasta_seqs)
        
        direct_matches = vcf_set.intersection(fasta_set)
        compatibility['potential_matches'] = len(direct_matches)
        
        if len(direct_matches) == 0:
            compatibility['warnings'].append(
                "No direct chromosome name matches between VCF and FASTA - "
                "this may be normal if chromosome naming conventions differ"
            )
        
        return compatibility
    
    def run_liftover(self):
        """Execute CrossMap liftover with quality validation"""
        logger.info("=== Step 2: CrossMap Liftover ===")
        
        logger.info(f"Running CrossMap liftover...")
        logger.info(f"  Chain: {self.chain_file}")
        logger.info(f"  Input VCF: {self.vcf_file}")
        logger.info(f"  Reference: {self.target_fasta}")
        logger.info(f"  Output: {self.lifted_vcf}")
        
        # Construct CrossMap command
        cmd = [
            'crossmap', 'vcf',
            str(self.chain_file),
            str(self.vcf_file),
            str(self.target_fasta),
            str(self.lifted_vcf)
        ]
        
        try:
            # Run CrossMap
            logger.info(f"Executing: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"CrossMap failed with return code {result.returncode}")
                logger.error(f"STDERR: {result.stderr}")
                logger.error(f"STDOUT: {result.stdout}")
                raise subprocess.CalledProcessError(result.returncode, cmd)
            
            logger.info("CrossMap completed successfully")
            if result.stdout:
                logger.info(f"CrossMap output: {result.stdout}")
            
            # Validate liftover results
            validation_result = self._validate_liftover_results()
            self._save_validation_report('liftover', validation_result)
            
            if not validation_result['valid']:
                raise ValueError(f"Liftover validation failed: {validation_result.get('error', 'Unknown error')}")
                
            logger.info(f"✓ CrossMap liftover completed successfully")
            logger.info(f"  Lifted variants: {validation_result['lifted_count']:,}")
            logger.info(f"  Unlifted variants: {validation_result['unlifted_count']:,}")
            logger.info(f"  Success rate: {validation_result['success_rate']:.1f}%")
            
            return validation_result
            
        except subprocess.TimeoutExpired:
            logger.error("CrossMap timed out after 1 hour")
            raise
        except FileNotFoundError:
            logger.error("CrossMap not found. Install with: conda install -c bioconda crossmap")
            raise
        except Exception as e:
            logger.error(f"Error running CrossMap: {e}")
            raise
    
    def _validate_liftover_results(self):
        """Validate liftover output files - FIXED VERSION"""
        validation = {
            'step': 'liftover_validation',
            'timestamp': datetime.now().isoformat(),
            'valid': False,
            'lifted_count': 0,
            'unlifted_count': 0,
            'success_rate': 0,
            'chromosome_distribution': Counter(),
            'quality_metrics': {}
        }
        
        try:
            # Count lifted variants
            if self.lifted_vcf.exists():
                with open(self.lifted_vcf, 'r') as f:
                    for line in f:
                        if not line.startswith('#'):
                            validation['lifted_count'] += 1
                            # Track chromosome distribution
                            chrom = line.split('\t')[0]
                            validation['chromosome_distribution'][chrom] += 1
            
            # Count unlifted variants - FIXED: Look for .vcf.unmap file
            unlifted_file = self.lifted_vcf.with_suffix('.vcf.unmap')
            if unlifted_file.exists():
                with open(unlifted_file, 'r') as f:
                    validation['unlifted_count'] = sum(1 for line in f if line.strip() and not line.startswith('#'))
            else:
                # Alternative: check if there's a .unmap file with different naming
                alternative_unlifted = self.output_dir / f"{self.output_prefix}_lifted.vcf.unmap"
                if alternative_unlifted.exists():
                    with open(alternative_unlifted, 'r') as f:
                        validation['unlifted_count'] = sum(1 for line in f if line.strip() and not line.startswith('#'))
            
            # Calculate success rate
            total_variants = validation['lifted_count'] + validation['unlifted_count']
            if total_variants > 0:
                validation['success_rate'] = (validation['lifted_count'] / total_variants) * 100
            
            # Quality assessment
            validation['quality_metrics'] = self._assess_liftover_quality(validation['success_rate'])
            
            validation['valid'] = validation['lifted_count'] > 0
            
            if not validation['valid']:
                validation['error'] = 'No variants were successfully lifted'
                
        except Exception as e:
            validation['error'] = str(e)
            
        return validation
    
    def _assess_liftover_quality(self, success_rate):
        """Assess overall liftover quality"""
        if success_rate >= 85:
            quality = 'Excellent'
            recommendation = 'Results are highly reliable for downstream analysis'
        elif success_rate >= 70:
            quality = 'Good'
            recommendation = 'Results are suitable for most analyses'
        elif success_rate >= 50:
            quality = 'Moderate'
            recommendation = 'Consider parameter optimization or manual review'
        else:
            quality = 'Poor'
            recommendation = 'Check alignment quality and consider alternative approaches'
            
        return {
            'overall_quality': quality,
            'recommendation': recommendation,
            'success_rate': success_rate
        }
    
    def finalize_output(self):
        """Compress and index final VCF with validation"""
        logger.info("=== Step 3: Output Finalization ===")
        
        validation = {
            'step': 'output_finalization',
            'timestamp': datetime.now().isoformat(),
            'valid': False,
            'final_files': {}
        }
        
        try:
            # Check if bgzip and tabix are available
            bgzip_available = subprocess.run(['which', 'bgzip'], capture_output=True).returncode == 0
            tabix_available = subprocess.run(['which', 'tabix'], capture_output=True).returncode == 0
            
            if bgzip_available and tabix_available:
                # Compress VCF
                compressed_vcf = self.lifted_vcf.with_suffix('.vcf.gz')
                subprocess.run(['bgzip', '-c', str(self.lifted_vcf)], 
                             stdout=open(compressed_vcf, 'wb'), check=True)
                
                # Index compressed VCF
                subprocess.run(['tabix', '-p', 'vcf', str(compressed_vcf)], check=True)
                
                validation['final_files']['compressed_vcf'] = {
                    'path': str(compressed_vcf),
                    'size_mb': round(compressed_vcf.stat().st_size / (1024*1024), 2),
                    'indexed': compressed_vcf.with_suffix('.vcf.gz.tbi').exists()
                }
                
                logger.info(f"✓ VCF compressed and indexed: {compressed_vcf}")
            else:
                logger.warning("bgzip/tabix not available - VCF not compressed")
                validation['final_files']['uncompressed_vcf'] = {
                    'path': str(self.lifted_vcf),
                    'size_mb': round(self.lifted_vcf.stat().st_size / (1024*1024), 2),
                    'indexed': False
                }
            
            # Validate other output files
            for file_type, file_path in [
                ('lifted_vcf', self.lifted_vcf),
                ('validation_dir', self.validation_dir)
            ]:
                if file_path.exists():
                    if file_path.is_file():
                        size_mb = round(file_path.stat().st_size / (1024*1024), 2)
                    else:
                        size_mb = sum(f.stat().st_size for f in file_path.rglob('*') if f.is_file()) / (1024*1024)
                        size_mb = round(size_mb, 2)
                    
                    validation['final_files'][file_type] = {
                        'path': str(file_path),
                        'size_mb': size_mb,
                        'exists': True
                    }
            
            validation['valid'] = len(validation['final_files']) > 0
            
        except Exception as e:
            validation['error'] = str(e)
            
        self._save_validation_report('output_finalization', validation)
        
        if validation['valid']:
            logger.info("✓ Output finalization completed")
        
        return validation
    
    def generate_final_report(self):
        """Generate comprehensive final report"""
        logger.info("=== Step 4: Final Report Generation ===")
        
        # Collect all validation results
        all_validations = {}
        for report_file in self.validation_dir.glob('*.json'):
            with open(report_file, 'r') as f:
                data = json.load(f)
                all_validations[data['step']] = data
        
        # Generate summary report
        summary = {
            'pipeline_info': {
                'timestamp': datetime.now().isoformat(),
                'pipeline_type': 'simplified_liftover',
                'input_files': {
                    'chain_file': str(self.chain_file),
                    'vcf_file': str(self.vcf_file),
                    'target_fasta': str(self.target_fasta)
                },
                'output_prefix': self.output_prefix,
                'output_directory': str(self.output_dir)
            },
            'step_results': all_validations,
            'overall_success': all(v.get('valid', False) for v in all_validations.values()),
            'recommendations': []
        }
        
        # Add recommendations based on results
        if 'liftover_validation' in all_validations:
            liftover_data = all_validations['liftover_validation']
            success_rate = liftover_data.get('success_rate', 0)
            
            if success_rate >= 85:
                summary['recommendations'].append("Excellent liftover quality - results are highly reliable")
            elif success_rate >= 70:
                summary['recommendations'].append("Good liftover quality - suitable for most analyses")
            elif success_rate >= 50:
                summary['recommendations'].append("Moderate success rate - consider manual review of key variants")
            else:
                summary['recommendations'].append("Low success rate - check alignment quality and parameters")
        
        # Save comprehensive report
        report_file = self.validation_dir / 'final_report.json'
        with open(report_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        
        # Generate human-readable summary
        summary_file = self.validation_dir / 'pipeline_summary.txt'
        with open(summary_file, 'w') as f:
            f.write("Simplified A. lyrata to A. arenosa Liftover Pipeline Summary\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Execution Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Overall Success: {'✓ PASSED' if summary['overall_success'] else '✗ FAILED'}\n\n")
            
            # Input files summary
            f.write("Input Files:\n")
            f.write("-" * 20 + "\n")
            for key, path in summary['pipeline_info']['input_files'].items():
                f.write(f"  {key}: {path}\n")
            f.write("\n")
            
            # Step-by-step results
            f.write("Pipeline Steps:\n")
            f.write("-" * 20 + "\n")
            step_names = {
                'input_validation': 'Input File Validation',
                'liftover': 'CrossMap Liftover',
                'output_finalization': 'Output Finalization'
            }
            
            for step_key, step_name in step_names.items():
                if step_key in all_validations:
                    result = all_validations[step_key]
                    status = '✓ PASSED' if result.get('valid', False) else '✗ FAILED'
                    f.write(f"  {step_name}: {status}\n")
                    
                    # Add specific metrics
                    if step_key == 'liftover':
                        f.write(f"    - Lifted variants: {result.get('lifted_count', 0):,}\n")
                        f.write(f"    - Success rate: {result.get('success_rate', 0):.1f}%\n")
                        
            f.write("\n")
            
            # Recommendations
            if summary['recommendations']:
                f.write("Recommendations:\n")
                f.write("-" * 20 + "\n")
                for rec in summary['recommendations']:
                    f.write(f"• {rec}\n")
                f.write("\n")
            
            # Output files
            f.write("Output Files:\n")
            f.write("-" * 20 + "\n")
            f.write(f"  Lifted VCF: {self.lifted_vcf}\n")
            if self.lifted_vcf.with_suffix('.vcf.gz').exists():
                f.write(f"  Compressed VCF: {self.lifted_vcf.with_suffix('.vcf.gz')}\n")
            if self.unlifted_bed.exists():
                f.write(f"  Unlifted variants: {self.unlifted_bed}\n")
            f.write(f"  Validation reports: {self.validation_dir}\n")
        
        logger.info(f"✓ Final reports generated:")
        logger.info(f"  Detailed: {report_file}")
        logger.info(f"  Summary: {summary_file}")
        
        return summary
    
    def _open_file(self, file_path):
        """Open file with automatic gzip detection"""
        if str(file_path).endswith('.gz'):
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'r')
    
    def _save_validation_report(self, step_name, validation_data):
        """Save validation report to JSON file"""
        report_file = self.validation_dir / f"{step_name}_validation.json"
        with open(report_file, 'w') as f:
            json.dump(validation_data, f, indent=2, default=str)
        
        logger.debug(f"Validation report saved: {report_file}")
    
    def run_pipeline(self):
        """Execute the simplified liftover pipeline"""
        logger.info("=" * 80)
        logger.info("SIMPLIFIED ARABIDOPSIS LIFTOVER PIPELINE - A. lyrata to A. arenosa")
        logger.info("=" * 80)
        
        start_time = datetime.now()
        
        try:
            # Step 1: Validate input files
            self.validate_input_files()
            
            # Step 2: Run liftover
            self.run_liftover()
            
            # Step 3: Finalize output
            self.finalize_output()
            
            # Step 4: Generate reports
            final_report = self.generate_final_report()
            
            end_time = datetime.now()
            duration = end_time - start_time
            
            logger.info("=" * 80)
            logger.info("PIPELINE COMPLETED SUCCESSFULLY!")
            logger.info(f"Total execution time: {duration}")
            logger.info(f"Output directory: {self.output_dir}")
            logger.info("=" * 80)
            
            return final_report
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            logger.error("Check validation reports for detailed error information")
            raise


def main():
    """Main function - runs simplified pipeline with configuration from top of script"""
    
    # Use configuration from top of script
    config = INPUT_CONFIG.copy()
    
    # Configure logging level
    if config.get('verbose', False):
        logging.getLogger().setLevel(logging.DEBUG)
        logger.info("Verbose logging enabled")
    
    # Print configuration
    logger.info("=" * 80)
    logger.info("SIMPLIFIED ARABIDOPSIS LIFTOVER PIPELINE CONFIGURATION")
    logger.info("=" * 80)
    logger.info("Input files:")
    logger.info(f"  Chain file: {config['chain_file']}")
    logger.info(f"  VCF file: {config['vcf_file']}")
    logger.info(f"  Target FASTA: {config['target_fasta']}")
    logger.info(f"Output configuration:")
    logger.info(f"  Directory: {config['output_dir']}")
    logger.info(f"  Prefix: {config['output_prefix']}")
    logger.info("=" * 80)
    
    # Validate required tools
    required_tools = ['crossmap']
    optional_tools = ['bgzip', 'tabix']
    missing_required = []
    missing_optional = []
    
    for tool in required_tools:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing_required.append(tool)
    
    for tool in optional_tools:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing_optional.append(tool)
    
    if missing_required:
        logger.error(f"Missing required tools: {', '.join(missing_required)}")
        logger.error("Install with: conda install -c bioconda crossmap")
        sys.exit(1)
    
    if missing_optional:
        logger.warning(f"Missing optional tools: {', '.join(missing_optional)}")
        logger.warning("Install with: conda install -c bioconda htslib")
        logger.warning("Output VCF will not be compressed/indexed")
    
    # Basic file existence check
    for key, file_path in [('chain_file', config['chain_file']), 
                          ('vcf_file', config['vcf_file']), 
                          ('target_fasta', config['target_fasta'])]:
        if not Path(file_path).exists():
            logger.error(f"File not found: {file_path}")
            logger.error(f"Please update the {key} path in the INPUT_CONFIG section at the top of this script")
            sys.exit(1)
    
    # Run pipeline
    try:
        pipeline = SimplifiedLiftoverPipeline(config)
        final_report = pipeline.run_pipeline()
        
        # Print summary to console
        print("\n" + "=" * 60)
        print("PIPELINE SUMMARY")
        print("=" * 60)
        
        if final_report['overall_success']:
            print("✓ Pipeline completed successfully!")
            
            # Extract key metrics
            if 'liftover_validation' in final_report['step_results']:
                liftover_data = final_report['step_results']['liftover_validation']
                print(f"✓ Lifted variants: {liftover_data.get('lifted_count', 0):,}")
                print(f"✓ Success rate: {liftover_data.get('success_rate', 0):.1f}%")
                
                # Quality assessment
                success_rate = liftover_data.get('success_rate', 0)
                if success_rate >= 85:
                    print("✓ Quality: Excellent - results are highly reliable")
                elif success_rate >= 70:
                    print("✓ Quality: Good - suitable for most analyses")
                elif success_rate >= 50:
                    print("⚠ Quality: Moderate - consider manual review")
                else:
                    print("⚠ Quality: Poor - check alignment parameters")
        else:
            print("✗ Pipeline failed - check logs for details")
            
        print(f"\nDetailed reports: {pipeline.validation_dir}")
        print(f"Output files: {pipeline.output_dir}")
        print("\nNext steps:")
        print(f"1. Review {pipeline.validation_dir}/pipeline_summary.txt")
        print(f"2. Use {pipeline.output_dir}/{config['output_prefix']}_lifted.vcf for analysis")
        print(f"3. Check unlifted variants in {pipeline.output_dir}/{config['output_prefix']}_unlifted.bed")
        
    except KeyboardInterrupt:
        logger.info("Pipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Pipeline failed with error: {e}")
        logger.error("Check the log file and validation reports for detailed error information")
        sys.exit(1)


if __name__ == "__main__":
    main()