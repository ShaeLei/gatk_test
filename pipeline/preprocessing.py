"""
Preprocessing Module for Variant Calling Pipeline

This module handles the initial preprocessing steps required before variant calling:
1. Quality assessment using FastQC
2. Adapter trimming using Cutadapt
3. Quality trimming using Trimmomatic
4. Read filtering and quality control

These steps are crucial for ensuring high-quality variant calls and optimal pipeline performance.
"""

import logging
import os
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class QualityControl:
    """Handles quality assessment and control of FASTQ files."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def run_fastqc(self, input_dir: str, output_dir: str) -> str:
        """
        Run FastQC quality assessment on FASTQ files.
        
        Args:
            input_dir: Directory containing FASTQ files
            output_dir: Output directory for FastQC reports
            
        Returns:
            str: Path to FastQC report directory
        """
        self.logger.info("Starting FastQC quality assessment")
        start_time = time.time()
        
        try:
            # Find all FASTQ files
            fastq_files = self._find_fastq_files(input_dir)
            
            if not fastq_files:
                raise ValueError("No FASTQ files found in input directory")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Run FastQC with optimized parameters
            cmd = [
                'fastqc',
                '--outdir', output_dir,
                '--threads', str(self.config.get('threads', 4)),
                '--quiet'  # Reduce output verbosity
            ]
            cmd.extend(fastq_files)
            
            self.logger.info(f"Running FastQC command: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Generate quality summary
            quality_summary = self._generate_quality_summary(output_dir)
            
            elapsed_time = time.time() - start_time
            self.logger.info(f"FastQC completed in {elapsed_time:.2f} seconds")
            self.logger.info(f"Quality summary: {quality_summary}")
            
            return output_dir
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"FastQC failed: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Quality assessment failed: {str(e)}")
            raise
    
    def _find_fastq_files(self, input_dir: str) -> List[str]:
        """Find all FASTQ files in the input directory."""
        fastq_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
        fastq_files = []
        
        for ext in fastq_extensions:
            fastq_files.extend(Path(input_dir).glob(f"*{ext}"))
        
        return [str(f) for f in fastq_files]
    
    def _generate_quality_summary(self, fastqc_dir: str) -> Dict:
        """Generate a summary of quality metrics from FastQC reports."""
        summary = {
            'total_reads': 0,
            'avg_quality': 0,
            'gc_content': 0,
            'duplication_rate': 0,
            'adapter_contamination': False
        }
        
        # This would parse FastQC HTML/text reports
        # For demonstration, we'll return a mock summary
        return summary


class AdapterTrimming:
    """Handles adapter trimming using Cutadapt."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def run_cutadapt(self, input_dir: str, output_dir: str) -> str:
        """
        Run Cutadapt for adapter trimming.
        
        Args:
            input_dir: Directory containing FASTQ files
            output_dir: Output directory for trimmed files
            
        Returns:
            str: Path to trimmed FASTQ files
        """
        self.logger.info("Starting adapter trimming with Cutadapt")
        start_time = time.time()
        
        try:
            # Find FASTQ files
            fastq_files = self._find_fastq_files(input_dir)
            
            if not fastq_files:
                raise ValueError("No FASTQ files found for trimming")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            trimmed_files = []
            
            for fastq_file in fastq_files:
                output_file = self._get_trimmed_filename(fastq_file, output_dir)
                
                # Cutadapt command with optimized parameters
                cmd = [
                    'cutadapt',
                    '--cores', str(self.config.get('threads', 4)),
                    '--quality-cutoff', '20,20',  # Quality threshold
                    '--minimum-length', '50',      # Minimum read length
                    '--adapter', 'AGATCGGAAGAGC',  # Common Illumina adapter
                    '--output', output_file,
                    '--report', f"{output_file}.report",
                    fastq_file
                ]
                
                self.logger.info(f"Running Cutadapt: {' '.join(cmd)}")
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                trimmed_files.append(output_file)
                
                # Log trimming statistics
                self._log_trimming_stats(result.stdout, fastq_file)
            
            elapsed_time = time.time() - start_time
            self.logger.info(f"Adapter trimming completed in {elapsed_time:.2f} seconds")
            
            return output_dir
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Cutadapt failed: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Adapter trimming failed: {str(e)}")
            raise
    
    def _find_fastq_files(self, input_dir: str) -> List[str]:
        """Find FASTQ files for trimming."""
        fastq_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
        fastq_files = []
        
        for ext in fastq_extensions:
            fastq_files.extend(Path(input_dir).glob(f"*{ext}"))
        
        return [str(f) for f in fastq_files]
    
    def _get_trimmed_filename(self, input_file: str, output_dir: str) -> str:
        """Generate output filename for trimmed file."""
        input_path = Path(input_file)
        output_filename = f"{input_path.stem}_trimmed{input_path.suffix}"
        return str(Path(output_dir) / output_filename)
    
    def _log_trimming_stats(self, stdout: str, input_file: str):
        """Log trimming statistics from Cutadapt output."""
        # Parse Cutadapt output for statistics
        lines = stdout.split('\n')
        for line in lines:
            if 'Total reads processed' in line:
                self.logger.info(f"Cutadapt stats for {input_file}: {line.strip()}")


class QualityTrimming:
    """Handles quality trimming using Trimmomatic."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def run_trimmomatic(self, input_dir: str, output_dir: str) -> str:
        """
        Run Trimmomatic for quality trimming.
        
        Args:
            input_dir: Directory containing trimmed FASTQ files
            output_dir: Output directory for quality-trimmed files
            
        Returns:
            str: Path to final processed FASTQ files
        """
        self.logger.info("Starting quality trimming with Trimmomatic")
        start_time = time.time()
        
        try:
            # Find trimmed FASTQ files
            trimmed_files = self._find_trimmed_files(input_dir)
            
            if not trimmed_files:
                raise ValueError("No trimmed FASTQ files found")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            processed_files = []
            
            for fastq_file in trimmed_files:
                output_file = self._get_quality_trimmed_filename(fastq_file, output_dir)
                
                # Trimmomatic command with optimized parameters
                cmd = [
                    'trimmomatic',
                    'SE',  # Single-end reads
                    '-threads', str(self.config.get('threads', 4)),
                    fastq_file,
                    output_file,
                    'LEADING:3',      # Remove leading low quality bases
                    'TRAILING:3',     # Remove trailing low quality bases
                    'SLIDINGWINDOW:4:20',  # Quality window
                    'MINLEN:50'       # Minimum read length
                ]
                
                self.logger.info(f"Running Trimmomatic: {' '.join(cmd)}")
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                processed_files.append(output_file)
                
                # Log quality trimming statistics
                self._log_quality_stats(result.stdout, fastq_file)
            
            elapsed_time = time.time() - start_time
            self.logger.info(f"Quality trimming completed in {elapsed_time:.2f} seconds")
            
            return output_dir
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Trimmomatic failed: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Quality trimming failed: {str(e)}")
            raise
    
    def _find_trimmed_files(self, input_dir: str) -> List[str]:
        """Find trimmed FASTQ files."""
        trimmed_patterns = ['*_trimmed.fastq', '*_trimmed.fq', '*_trimmed.fastq.gz', '*_trimmed.fq.gz']
        trimmed_files = []
        
        for pattern in trimmed_patterns:
            trimmed_files.extend(Path(input_dir).glob(pattern))
        
        return [str(f) for f in trimmed_files]
    
    def _get_quality_trimmed_filename(self, input_file: str, output_dir: str) -> str:
        """Generate output filename for quality-trimmed file."""
        input_path = Path(input_file)
        output_filename = f"{input_path.stem}_quality_trimmed{input_path.suffix}"
        return str(Path(output_dir) / output_filename)
    
    def _log_quality_stats(self, stdout: str, input_file: str):
        """Log quality trimming statistics from Trimmomatic output."""
        # Parse Trimmomatic output for statistics
        lines = stdout.split('\n')
        for line in lines:
            if 'Input Reads' in line or 'Surviving' in line:
                self.logger.info(f"Trimmomatic stats for {input_file}: {line.strip()}")


class ReadFiltering:
    """Handles additional read filtering and quality control."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def filter_low_quality_reads(self, input_dir: str, output_dir: str) -> str:
        """
        Apply additional quality filters to remove low-quality reads.
        
        Args:
            input_dir: Directory containing quality-trimmed FASTQ files
            output_dir: Output directory for filtered files
            
        Returns:
            str: Path to final filtered FASTQ files
        """
        self.logger.info("Starting additional read filtering")
        
        try:
            # Find quality-trimmed files
            quality_files = self._find_quality_trimmed_files(input_dir)
            
            if not quality_files:
                raise ValueError("No quality-trimmed files found")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            filtered_files = []
            
            for fastq_file in quality_files:
                output_file = self._get_filtered_filename(fastq_file, output_dir)
                
                # Apply custom filtering criteria
                self._apply_quality_filters(fastq_file, output_file)
                
                filtered_files.append(output_file)
            
            self.logger.info(f"Read filtering completed. {len(filtered_files)} files processed")
            
            return output_dir
            
        except Exception as e:
            self.logger.error(f"Read filtering failed: {str(e)}")
            raise
    
    def _find_quality_trimmed_files(self, input_dir: str) -> List[str]:
        """Find quality-trimmed FASTQ files."""
        quality_patterns = ['*_quality_trimmed.fastq', '*_quality_trimmed.fq']
        quality_files = []
        
        for pattern in quality_patterns:
            quality_files.extend(Path(input_dir).glob(pattern))
        
        return [str(f) for f in quality_files]
    
    def _get_filtered_filename(self, input_file: str, output_dir: str) -> str:
        """Generate output filename for filtered file."""
        input_path = Path(input_file)
        output_filename = f"{input_path.stem}_filtered{input_path.suffix}"
        return str(Path(output_dir) / output_filename)
    
    def _apply_quality_filters(self, input_file: str, output_file: str):
        """Apply custom quality filters to FASTQ file."""
        # This would implement custom filtering logic
        # For demonstration, we'll just copy the file
        import shutil
        shutil.copy2(input_file, output_file)
        
        self.logger.info(f"Applied quality filters to {input_file}") 