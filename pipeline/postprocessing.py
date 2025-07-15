"""
Post-processing Module for Variant Calling Pipeline

This module handles GATK post-processing steps that improve the quality
of aligned reads before variant calling:

1. MarkDuplicates: Remove PCR duplicates
2. BaseRecalibrator: Recalibrate base quality scores
3. ApplyBQSR: Apply recalibrated quality scores

These steps are crucial for accurate variant calling following GATK best practices.
"""

import logging
import os
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional


class Deduplication:
    """Handles PCR duplicate removal using GATK MarkDuplicates."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.gatk_path = config.get('gatk_path', 'gatk')
    
    def run_mark_duplicates(self, bam_file: str, output_dir: str) -> str:
        """
        Run GATK MarkDuplicates to remove PCR duplicates.
        
        Args:
            bam_file: Path to input BAM file
            output_dir: Output directory for deduplicated BAM
            
        Returns:
            str: Path to deduplicated BAM file
        """
        self.logger.info("Starting PCR duplicate removal with GATK MarkDuplicates")
        start_time = time.time()
        
        try:
            # Validate input
            if not os.path.exists(bam_file):
                raise ValueError(f"Input BAM file not found: {bam_file}")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Generate output filenames
            dedup_bam = os.path.join(output_dir, "deduplicated.bam")
            metrics_file = os.path.join(output_dir, "duplicate_metrics.txt")
            
            # Build GATK MarkDuplicates command with optimized parameters
            cmd = [
                self.gatk_path, 'MarkDuplicates',
                '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
                '--input', bam_file,
                '--output', dedup_bam,
                '--metrics', metrics_file,
                '--create-index', 'true',
                '--validation-stringency', 'SILENT'
            ]
            
            self.logger.info(f"Running GATK MarkDuplicates: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Log deduplication statistics
            self._log_deduplication_stats(metrics_file)
            
            elapsed_time = time.time() - start_time
            self.logger.info(f"MarkDuplicates completed in {elapsed_time:.2f} seconds")
            
            return dedup_bam
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"GATK MarkDuplicates failed: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Deduplication failed: {str(e)}")
            raise
    
    def _log_deduplication_stats(self, metrics_file: str):
        """Log deduplication statistics from metrics file."""
        try:
            with open(metrics_file, 'r') as f:
                lines = f.readlines()
            
            # Parse key metrics
            for line in lines:
                if line.startswith('LIBRARY') or line.startswith('UNPAIRED_READ_DUPLICATES'):
                    self.logger.info(f"Deduplication stats: {line.strip()}")
                    
        except Exception as e:
            self.logger.warning(f"Could not read deduplication metrics: {str(e)}")


class BaseRecalibration:
    """Handles base quality score recalibration using GATK BaseRecalibrator."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.gatk_path = config.get('gatk_path', 'gatk')
    
    def run_base_recalibrator(self, bam_file: str, reference: str, output_dir: str) -> str:
        """
        Run GATK BaseRecalibrator to recalibrate base quality scores.
        
        Args:
            bam_file: Path to deduplicated BAM file
            reference: Path to reference genome
            output_dir: Output directory for recalibrated BAM
            
        Returns:
            str: Path to recalibrated BAM file
        """
        self.logger.info("Starting base quality score recalibration")
        start_time = time.time()
        
        try:
            # Validate inputs
            if not os.path.exists(bam_file):
                raise ValueError(f"Input BAM file not found: {bam_file}")
            if not os.path.exists(reference):
                raise ValueError(f"Reference genome not found: {reference}")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Step 1: Generate recalibration table
            recal_table = self._generate_recalibration_table(bam_file, reference, output_dir)
            
            # Step 2: Apply recalibration
            recal_bam = self._apply_recalibration(bam_file, recal_table, output_dir)
            
            elapsed_time = time.time() - start_time
            self.logger.info(f"Base recalibration completed in {elapsed_time:.2f} seconds")
            
            return recal_bam
            
        except Exception as e:
            self.logger.error(f"Base recalibration failed: {str(e)}")
            raise
    
    def _generate_recalibration_table(self, bam_file: str, reference: str, output_dir: str) -> str:
        """Generate recalibration table using GATK BaseRecalibrator."""
        recal_table = os.path.join(output_dir, "recalibration.table")
        
        # Known sites for recalibration (these would be provided)
        known_sites = self._get_known_sites()
        
        cmd = [
            self.gatk_path, 'BaseRecalibrator',
            '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
            '--input', bam_file,
            '--reference', reference,
            '--output', recal_table
        ]
        
        # Add known sites if available
        for site in known_sites:
            cmd.extend(['--known-sites', site])
        
        self.logger.info(f"Generating recalibration table: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        return recal_table
    
    def _apply_recalibration(self, bam_file: str, recal_table: str, output_dir: str) -> str:
        """Apply recalibration using GATK ApplyBQSR."""
        recal_bam = os.path.join(output_dir, "recalibrated.bam")
        
        cmd = [
            self.gatk_path, 'ApplyBQSR',
            '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
            '--input', bam_file,
            '--bqsr-recal-file', recal_table,
            '--output', recal_bam,
            '--create-index', 'true'
        ]
        
        self.logger.info(f"Applying recalibration: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        return recal_bam
    
    def _get_known_sites(self) -> List[str]:
        """Get known variant sites for recalibration."""
        # In a real implementation, these would be paths to known variant databases
        # such as dbSNP, 1000 Genomes, etc.
        known_sites = []
        
        # Example known sites (these would be actual file paths)
        # known_sites = [
        #     "/path/to/dbsnp_138.hg38.vcf",
        #     "/path/to/1000G_phase1.snps.high_confidence.hg38.vcf"
        # ]
        
        return known_sites


class QualityImprovement:
    """Handles additional quality improvement steps."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.gatk_path = config.get('gatk_path', 'gatk')
    
    def run_quality_improvement(self, bam_file: str, reference: str, output_dir: str) -> str:
        """
        Run additional quality improvement steps.
        
        Args:
            bam_file: Path to recalibrated BAM file
            reference: Path to reference genome
            output_dir: Output directory for improved BAM
            
        Returns:
            str: Path to final processed BAM file
        """
        self.logger.info("Starting additional quality improvement steps")
        
        try:
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Step 1: Fix mate information (for paired-end reads)
            fixed_bam = self._fix_mate_information(bam_file, output_dir)
            
            # Step 2: Validate BAM file
            self._validate_bam(fixed_bam, reference)
            
            # Step 3: Generate final statistics
            self._generate_final_stats(fixed_bam, output_dir)
            
            return fixed_bam
            
        except Exception as e:
            self.logger.error(f"Quality improvement failed: {str(e)}")
            raise
    
    def _fix_mate_information(self, bam_file: str, output_dir: str) -> str:
        """Fix mate information for paired-end reads."""
        fixed_bam = os.path.join(output_dir, "fixed_mates.bam")
        
        cmd = [
            self.gatk_path, 'FixMateInformation',
            '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
            '--input', bam_file,
            '--output', fixed_bam,
            '--create-index', 'true'
        ]
        
        self.logger.info(f"Fixing mate information: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        return fixed_bam
    
    def _validate_bam(self, bam_file: str, reference: str):
        """Validate BAM file using GATK ValidateSamFile."""
        cmd = [
            self.gatk_path, 'ValidateSamFile',
            '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
            '--input', bam_file,
            '--reference', reference,
            '--mode', 'VERBOSE'
        ]
        
        self.logger.info(f"Validating BAM file: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Check for validation errors
        if 'ERROR' in result.stdout or 'ERROR' in result.stderr:
            self.logger.warning("BAM validation found errors")
        else:
            self.logger.info("BAM validation passed")
    
    def _generate_final_stats(self, bam_file: str, output_dir: str):
        """Generate final BAM statistics."""
        stats_file = os.path.join(output_dir, "final_bam_stats.txt")
        
        cmd = [
            self.gatk_path, 'FlagStat',
            '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
            '--input', bam_file,
            '--output', stats_file
        ]
        
        self.logger.info(f"Generating final BAM statistics: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Log key statistics
        self._log_final_stats(stats_file)
    
    def _log_final_stats(self, stats_file: str):
        """Log final BAM statistics."""
        try:
            with open(stats_file, 'r') as f:
                stats = f.read()
            
            lines = stats.split('\n')
            for line in lines:
                if 'total' in line or 'mapped' in line or 'paired' in line:
                    self.logger.info(f"Final BAM stats: {line.strip()}")
                    
        except Exception as e:
            self.logger.warning(f"Could not read final stats: {str(e)}")


class PostProcessingOptimizer:
    """Handles performance optimization for post-processing steps."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def optimize_parameters(self, bam_file: str) -> Dict:
        """
        Optimize parameters based on BAM file characteristics.
        
        Args:
            bam_file: Path to BAM file for analysis
            
        Returns:
            Dict: Optimized parameters
        """
        self.logger.info("Optimizing post-processing parameters")
        
        try:
            # Analyze BAM file characteristics
            file_size = os.path.getsize(bam_file) / (1024 * 1024 * 1024)  # GB
            read_count = self._estimate_read_count(bam_file)
            
            # Optimize memory allocation
            optimal_memory = self._calculate_optimal_memory(file_size, read_count)
            
            # Optimize thread allocation
            optimal_threads = self._calculate_optimal_threads()
            
            optimized_params = {
                'memory': optimal_memory,
                'threads': optimal_threads,
                'chunk_size': self._calculate_chunk_size(file_size)
            }
            
            self.logger.info(f"Optimized parameters: {optimized_params}")
            
            return optimized_params
            
        except Exception as e:
            self.logger.error(f"Parameter optimization failed: {str(e)}")
            # Return default parameters
            return {
                'memory': self.config.get('memory', 8),
                'threads': self.config.get('threads', 4),
                'chunk_size': 1000000
            }
    
    def _estimate_read_count(self, bam_file: str) -> int:
        """Estimate number of reads in BAM file."""
        # This would use samtools to count reads
        # For demonstration, returning a mock value
        return 10000000  # 10M reads
    
    def _calculate_optimal_memory(self, file_size_gb: float, read_count: int) -> int:
        """Calculate optimal memory allocation based on file size and read count."""
        # Simple heuristic: 1GB per 5GB of BAM file, minimum 4GB
        base_memory = max(4, int(file_size_gb / 5))
        
        # Cap at available system memory (assuming 32GB system)
        max_memory = 24  # Leave some memory for system
        
        return min(base_memory, max_memory)
    
    def _calculate_optimal_threads(self) -> int:
        """Calculate optimal thread allocation."""
        # Use 75% of available CPU cores
        import multiprocessing
        available_cores = multiprocessing.cpu_count()
        optimal_threads = max(1, int(available_cores * 0.75))
        
        return optimal_threads
    
    def _calculate_chunk_size(self, file_size_gb: float) -> int:
        """Calculate optimal chunk size for processing."""
        # Larger files benefit from larger chunks
        if file_size_gb > 50:
            return 5000000  # 5M reads per chunk
        elif file_size_gb > 10:
            return 2000000  # 2M reads per chunk
        else:
            return 1000000  # 1M reads per chunk 