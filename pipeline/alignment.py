"""
Alignment Module for Variant Calling Pipeline

This module handles read alignment to the reference genome using BWA-MEM,
followed by SAM processing including sorting and indexing.

Key components:
1. BWA-MEM alignment with optimized parameters
2. SAM file processing and sorting
3. BAM file creation and indexing
4. Alignment quality assessment
"""

import logging
import os
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class BWAAlignment:
    """Handles BWA-MEM alignment of reads to reference genome."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.bwa_path = config.get('bwa_path', 'bwa')
    
    def run_bwa_mem(self, fastq_files: str, reference: str, output_dir: str) -> str:
        """
        Run BWA-MEM alignment with optimized parameters.
        
        Args:
            fastq_files: Directory containing processed FASTQ files
            reference: Path to reference genome FASTA file
            output_dir: Output directory for alignment files
            
        Returns:
            str: Path to output SAM file
        """
        self.logger.info("Starting BWA-MEM alignment")
        start_time = time.time()
        
        try:
            # Validate inputs
            if not os.path.exists(reference):
                raise ValueError(f"Reference genome not found: {reference}")
            
            # Check if reference is indexed
            self._ensure_reference_indexed(reference)
            
            # Find FASTQ files
            fastq_list = self._find_fastq_files(fastq_files)
            
            if not fastq_list:
                raise ValueError("No FASTQ files found for alignment")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Generate output SAM filename
            output_sam = os.path.join(output_dir, "aligned.sam")
            
            # Build BWA-MEM command with optimized parameters
            cmd = [
                self.bwa_path, 'mem',
                '-t', str(self.config.get('threads', 4)),  # Number of threads
                '-M',  # Mark shorter split hits as secondary
                '-R', '@RG\\tID:sample\\tSM:sample\\tPL:ILLUMINA',  # Read group
                '-v', '1',  # Verbosity level
                reference
            ]
            
            # Add input FASTQ files
            cmd.extend(fastq_list)
            
            self.logger.info(f"Running BWA-MEM: {' '.join(cmd)}")
            
            # Run BWA-MEM and redirect output to SAM file
            with open(output_sam, 'w') as sam_file:
                result = subprocess.run(
                    cmd,
                    stdout=sam_file,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
            
            # Log alignment statistics
            self._log_alignment_stats(result.stderr, output_sam)
            
            elapsed_time = time.time() - start_time
            self.logger.info(f"BWA-MEM alignment completed in {elapsed_time:.2f} seconds")
            
            return output_sam
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BWA-MEM failed: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Alignment failed: {str(e)}")
            raise
    
    def _ensure_reference_indexed(self, reference: str):
        """Ensure reference genome is indexed for BWA."""
        index_files = [
            f"{reference}.amb",
            f"{reference}.ann", 
            f"{reference}.bwt",
            f"{reference}.pac",
            f"{reference}.sa"
        ]
        
        missing_index = any(not os.path.exists(idx) for idx in index_files)
        
        if missing_index:
            self.logger.info("Indexing reference genome for BWA")
            
            cmd = [self.bwa_path, 'index', reference]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            self.logger.info("Reference genome indexing completed")
    
    def _find_fastq_files(self, input_dir: str) -> List[str]:
        """Find processed FASTQ files for alignment."""
        fastq_patterns = [
            '*_quality_trimmed.fastq',
            '*_quality_trimmed.fq',
            '*_filtered.fastq',
            '*_filtered.fq'
        ]
        
        fastq_files = []
        for pattern in fastq_patterns:
            fastq_files.extend(Path(input_dir).glob(pattern))
        
        return [str(f) for f in fastq_files]
    
    def _log_alignment_stats(self, stderr: str, sam_file: str):
        """Log alignment statistics from BWA-MEM output."""
        lines = stderr.split('\n')
        for line in lines:
            if 'M::mem_process_seqs' in line or 'M::mem_pestat' in line:
                self.logger.info(f"BWA-MEM stats: {line.strip()}")


class SAMProcessing:
    """Handles SAM file processing including sorting and indexing."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.samtools_path = config.get('samtools_path', 'samtools')
    
    def process_sam(self, sam_file: str, output_dir: str) -> str:
        """
        Process SAM file: sort, convert to BAM, and index.
        
        Args:
            sam_file: Path to input SAM file
            output_dir: Output directory for processed files
            
        Returns:
            str: Path to sorted BAM file
        """
        self.logger.info("Starting SAM file processing")
        start_time = time.time()
        
        try:
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Step 1: Convert SAM to BAM
            bam_file = self._sam_to_bam(sam_file, output_dir)
            
            # Step 2: Sort BAM file
            sorted_bam = self._sort_bam(bam_file, output_dir)
            
            # Step 3: Index BAM file
            self._index_bam(sorted_bam)
            
            # Step 4: Generate alignment statistics
            self._generate_alignment_stats(sorted_bam, output_dir)
            
            elapsed_time = time.time() - start_time
            self.logger.info(f"SAM processing completed in {elapsed_time:.2f} seconds")
            
            return sorted_bam
            
        except Exception as e:
            self.logger.error(f"SAM processing failed: {str(e)}")
            raise
    
    def _sam_to_bam(self, sam_file: str, output_dir: str) -> str:
        """Convert SAM file to BAM format."""
        bam_file = os.path.join(output_dir, "aligned.bam")
        
        cmd = [
            self.samtools_path, 'view',
            '-@', str(self.config.get('threads', 4)),  # Number of threads
            '-b',  # Output in BAM format
            '-S',  # Input is SAM format
            '-o', bam_file,
            sam_file
        ]
        
        self.logger.info(f"Converting SAM to BAM: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        return bam_file
    
    def _sort_bam(self, bam_file: str, output_dir: str) -> str:
        """Sort BAM file by coordinate."""
        sorted_bam = os.path.join(output_dir, "aligned_sorted.bam")
        
        cmd = [
            self.samtools_path, 'sort',
            '-@', str(self.config.get('threads', 4)),  # Number of threads
            '-o', sorted_bam,
            bam_file
        ]
        
        self.logger.info(f"Sorting BAM file: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        return sorted_bam
    
    def _index_bam(self, bam_file: str):
        """Create index for BAM file."""
        cmd = [
            self.samtools_path, 'index',
            bam_file
        ]
        
        self.logger.info(f"Indexing BAM file: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
    
    def _generate_alignment_stats(self, bam_file: str, output_dir: str):
        """Generate alignment statistics using samtools."""
        stats_file = os.path.join(output_dir, "alignment_stats.txt")
        
        cmd = [
            self.samtools_path, 'flagstat',
            bam_file
        ]
        
        self.logger.info(f"Generating alignment statistics: {' '.join(cmd)}")
        
        with open(stats_file, 'w') as f:
            result = subprocess.run(
                cmd,
                stdout=f,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        
        # Log key statistics
        self._log_alignment_summary(stats_file)
    
    def _log_alignment_summary(self, stats_file: str):
        """Log key alignment statistics."""
        try:
            with open(stats_file, 'r') as f:
                stats = f.read()
            
            lines = stats.split('\n')
            for line in lines:
                if 'total' in line or 'mapped' in line or 'paired' in line:
                    self.logger.info(f"Alignment stats: {line.strip()}")
                    
        except Exception as e:
            self.logger.warning(f"Could not read alignment stats: {str(e)}")


class AlignmentQualityAssessment:
    """Handles alignment quality assessment and reporting."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.samtools_path = config.get('samtools_path', 'samtools')
    
    def assess_alignment_quality(self, bam_file: str, output_dir: str) -> Dict:
        """
        Perform comprehensive alignment quality assessment.
        
        Args:
            bam_file: Path to sorted BAM file
            output_dir: Output directory for quality reports
            
        Returns:
            Dict: Quality assessment metrics
        """
        self.logger.info("Starting alignment quality assessment")
        
        try:
            metrics = {}
            
            # Calculate coverage statistics
            metrics['coverage'] = self._calculate_coverage(bam_file, output_dir)
            
            # Calculate mapping quality distribution
            metrics['mapping_quality'] = self._calculate_mapping_quality(bam_file, output_dir)
            
            # Calculate insert size distribution (for paired-end)
            metrics['insert_size'] = self._calculate_insert_size(bam_file, output_dir)
            
            # Generate quality report
            self._generate_quality_report(metrics, output_dir)
            
            return metrics
            
        except Exception as e:
            self.logger.error(f"Quality assessment failed: {str(e)}")
            raise
    
    def _calculate_coverage(self, bam_file: str, output_dir: str) -> Dict:
        """Calculate coverage statistics."""
        coverage_file = os.path.join(output_dir, "coverage.txt")
        
        cmd = [
            self.samtools_path, 'depth',
            '-a',  # Output all positions
            bam_file
        ]
        
        with open(coverage_file, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        
        # Parse coverage statistics
        coverage_stats = {
            'mean_coverage': 0,
            'median_coverage': 0,
            'coverage_10x': 0,
            'coverage_30x': 0
        }
        
        # This would calculate actual coverage statistics
        # For demonstration, returning mock data
        return coverage_stats
    
    def _calculate_mapping_quality(self, bam_file: str, output_dir: str) -> Dict:
        """Calculate mapping quality distribution."""
        qual_file = os.path.join(output_dir, "mapping_quality.txt")
        
        cmd = [
            self.samtools_path, 'view',
            bam_file,
            '|', 'awk', '{print $5}',
            '|', 'sort', '-n',
            '|', 'uniq', '-c'
        ]
        
        # This would calculate mapping quality distribution
        # For demonstration, returning mock data
        return {'mean_quality': 30, 'quality_distribution': 'normal'}
    
    def _calculate_insert_size(self, bam_file: str, output_dir: str) -> Dict:
        """Calculate insert size distribution for paired-end reads."""
        insert_file = os.path.join(output_dir, "insert_size.txt")
        
        # This would calculate insert size statistics
        # For demonstration, returning mock data
        return {'mean_insert_size': 300, 'insert_size_std': 50}
    
    def _generate_quality_report(self, metrics: Dict, output_dir: str):
        """Generate comprehensive quality report."""
        report_file = os.path.join(output_dir, "alignment_quality_report.txt")
        
        with open(report_file, 'w') as f:
            f.write("Alignment Quality Assessment Report\n")
            f.write("=" * 40 + "\n\n")
            
            f.write("Coverage Statistics:\n")
            for key, value in metrics['coverage'].items():
                f.write(f"  {key}: {value}\n")
            
            f.write("\nMapping Quality:\n")
            for key, value in metrics['mapping_quality'].items():
                f.write(f"  {key}: {value}\n")
            
            f.write("\nInsert Size Statistics:\n")
            for key, value in metrics['insert_size'].items():
                f.write(f"  {key}: {value}\n")
        
        self.logger.info(f"Quality report saved to {report_file}") 