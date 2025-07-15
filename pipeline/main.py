#!/usr/bin/env python3
"""
Variant Calling Pipeline - Main Orchestration Script

This script implements a complete variant calling pipeline following GATK best practices.
It demonstrates the workflow from raw FASTQ files to final variant calls with
performance optimization and accuracy validation.
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

# Import pipeline modules
from preprocessing import QualityControl, AdapterTrimming, QualityTrimming
from alignment import BWAAlignment, SAMProcessing
from postprocessing import Deduplication, BaseRecalibration
from variant_calling import HaplotypeCaller, VariantFiltering
from validation import ValidationMetrics, PerformanceMonitor


class VariantCallingPipeline:
    """
    Main pipeline class that orchestrates the complete variant calling workflow.
    
    This class demonstrates:
    1. Proper workflow design and modularity
    2. Performance monitoring and optimization
    3. Error handling and recovery
    4. Validation and quality assessment
    """
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = self._setup_logging()
        self.performance_monitor = PerformanceMonitor()
        
        # Initialize pipeline components
        self.qc = QualityControl(config)
        self.adapter_trim = AdapterTrimming(config)
        self.quality_trim = QualityTrimming(config)
        self.alignment = BWAAlignment(config)
        self.sam_processing = SAMProcessing(config)
        self.deduplication = Deduplication(config)
        self.recalibration = BaseRecalibration(config)
        self.variant_caller = HaplotypeCaller(config)
        self.variant_filter = VariantFiltering(config)
        self.validation = ValidationMetrics(config)
    
    def _setup_logging(self) -> logging.Logger:
        """Setup comprehensive logging for the pipeline."""
        logger = logging.getLogger('VariantCallingPipeline')
        logger.setLevel(logging.INFO)
        
        # Console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        
        # File handler
        fh = logging.FileHandler('pipeline.log')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        
        return logger
    
    def run_pipeline(self, input_dir: str, reference: str, output_dir: str) -> bool:
        """
        Execute the complete variant calling pipeline.
        
        Args:
            input_dir: Directory containing FASTQ files
            reference: Path to reference genome
            output_dir: Output directory for results
            
        Returns:
            bool: True if pipeline completed successfully
        """
        try:
            self.logger.info("Starting variant calling pipeline")
            self.performance_monitor.start_pipeline()
            
            # Create output directory structure
            self._setup_output_dirs(output_dir)
            
            # Step 1: Preprocessing
            self.logger.info("Step 1: Quality Control and Preprocessing")
            processed_fastq = self._run_preprocessing(input_dir, output_dir)
            
            # Step 2: Alignment
            self.logger.info("Step 2: Read Alignment")
            aligned_bam = self._run_alignment(processed_fastq, reference, output_dir)
            
            # Step 3: Post-processing
            self.logger.info("Step 3: Post-processing and Quality Improvement")
            processed_bam = self._run_postprocessing(aligned_bam, reference, output_dir)
            
            # Step 4: Variant Calling
            self.logger.info("Step 4: Variant Calling")
            raw_vcf = self._run_variant_calling(processed_bam, output_dir)
            
            # Step 5: Variant Filtering
            self.logger.info("Step 5: Variant Filtering")
            filtered_vcf = self._run_variant_filtering(raw_vcf, output_dir)
            
            # Step 6: Validation and Quality Assessment
            self.logger.info("Step 6: Validation and Quality Assessment")
            self._run_validation(filtered_vcf, output_dir)
            
            # Performance summary
            self.performance_monitor.end_pipeline()
            self._generate_performance_report(output_dir)
            
            self.logger.info("Pipeline completed successfully!")
            return True
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            self.performance_monitor.record_error(str(e))
            return False
    
    def _setup_output_dirs(self, output_dir: str):
        """Create organized output directory structure."""
        dirs = [
            'preprocessing',
            'alignment', 
            'postprocessing',
            'variants',
            'validation',
            'reports'
        ]
        
        for dir_name in dirs:
            Path(f"{output_dir}/{dir_name}").mkdir(parents=True, exist_ok=True)
    
    def _run_preprocessing(self, input_dir: str, output_dir: str) -> str:
        """Execute preprocessing steps with performance monitoring."""
        self.performance_monitor.start_step("preprocessing")
        
        try:
            # Quality assessment
            self.logger.info("Running quality assessment...")
            qc_report = self.qc.run_fastqc(input_dir, f"{output_dir}/preprocessing")
            
            # Adapter trimming
            self.logger.info("Running adapter trimming...")
            trimmed_fastq = self.adapter_trim.run_cutadapt(
                input_dir, f"{output_dir}/preprocessing"
            )
            
            # Quality trimming
            self.logger.info("Running quality trimming...")
            final_fastq = self.quality_trim.run_trimmomatic(
                trimmed_fastq, f"{output_dir}/preprocessing"
            )
            
            self.performance_monitor.end_step("preprocessing")
            return final_fastq
            
        except Exception as e:
            self.logger.error(f"Preprocessing failed: {str(e)}")
            raise
    
    def _run_alignment(self, fastq_files: str, reference: str, output_dir: str) -> str:
        """Execute alignment steps with performance monitoring."""
        self.performance_monitor.start_step("alignment")
        
        try:
            # BWA alignment
            self.logger.info("Running BWA alignment...")
            sam_file = self.alignment.run_bwa_mem(
                fastq_files, reference, f"{output_dir}/alignment"
            )
            
            # SAM processing (sorting, indexing)
            self.logger.info("Processing SAM files...")
            bam_file = self.sam_processing.process_sam(
                sam_file, f"{output_dir}/alignment"
            )
            
            self.performance_monitor.end_step("alignment")
            return bam_file
            
        except Exception as e:
            self.logger.error(f"Alignment failed: {str(e)}")
            raise
    
    def _run_postprocessing(self, bam_file: str, reference: str, output_dir: str) -> str:
        """Execute post-processing steps with performance monitoring."""
        self.performance_monitor.start_step("postprocessing")
        
        try:
            # Deduplication
            self.logger.info("Running deduplication...")
            dedup_bam = self.deduplication.run_mark_duplicates(
                bam_file, f"{output_dir}/postprocessing"
            )
            
            # Base quality score recalibration
            self.logger.info("Running base quality score recalibration...")
            recal_bam = self.recalibration.run_base_recalibrator(
                dedup_bam, reference, f"{output_dir}/postprocessing"
            )
            
            self.performance_monitor.end_step("postprocessing")
            return recal_bam
            
        except Exception as e:
            self.logger.error(f"Post-processing failed: {str(e)}")
            raise
    
    def _run_variant_calling(self, bam_file: str, output_dir: str) -> str:
        """Execute variant calling with performance monitoring."""
        self.performance_monitor.start_step("variant_calling")
        
        try:
            # HaplotypeCaller
            self.logger.info("Running HaplotypeCaller...")
            raw_vcf = self.variant_caller.run_haplotype_caller(
                bam_file, f"{output_dir}/variants"
            )
            
            self.performance_monitor.end_step("variant_calling")
            return raw_vcf
            
        except Exception as e:
            self.logger.error(f"Variant calling failed: {str(e)}")
            raise
    
    def _run_variant_filtering(self, vcf_file: str, output_dir: str) -> str:
        """Execute variant filtering with performance monitoring."""
        self.performance_monitor.start_step("variant_filtering")
        
        try:
            # Variant filtering
            self.logger.info("Running variant filtering...")
            filtered_vcf = self.variant_filter.run_variant_filtration(
                vcf_file, f"{output_dir}/variants"
            )
            
            self.performance_monitor.end_step("variant_filtering")
            return filtered_vcf
            
        except Exception as e:
            self.logger.error(f"Variant filtering failed: {str(e)}")
            raise
    
    def _run_validation(self, vcf_file: str, output_dir: str):
        """Execute validation and quality assessment."""
        self.performance_monitor.start_step("validation")
        
        try:
            # Calculate validation metrics
            self.logger.info("Calculating validation metrics...")
            metrics = self.validation.calculate_metrics(vcf_file)
            
            # Generate validation report
            self.logger.info("Generating validation report...")
            self.validation.generate_report(metrics, f"{output_dir}/validation")
            
            self.performance_monitor.end_step("validation")
            
        except Exception as e:
            self.logger.error(f"Validation failed: {str(e)}")
            raise
    
    def _generate_performance_report(self, output_dir: str):
        """Generate comprehensive performance report."""
        report = self.performance_monitor.generate_report()
        
        with open(f"{output_dir}/reports/performance_report.txt", 'w') as f:
            f.write("Variant Calling Pipeline Performance Report\n")
            f.write("=" * 50 + "\n\n")
            f.write(report)
        
        self.logger.info(f"Performance report saved to {output_dir}/reports/performance_report.txt")


def main():
    """Main entry point for the variant calling pipeline."""
    parser = argparse.ArgumentParser(
        description="Variant Calling Pipeline using GATK Best Practices"
    )
    parser.add_argument(
        "--input", required=True,
        help="Directory containing FASTQ files"
    )
    parser.add_argument(
        "--reference", required=True,
        help="Path to reference genome FASTA file"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--config", default="config/pipeline_config.yaml",
        help="Pipeline configuration file"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Number of threads for parallel processing"
    )
    parser.add_argument(
        "--memory", type=int, default=8,
        help="Memory allocation in GB"
    )
    
    args = parser.parse_args()
    
    # Load configuration
    config = {
        'threads': args.threads,
        'memory': args.memory,
        'gatk_path': os.getenv('GATK_PATH', '/usr/local/bin/gatk'),
        'bwa_path': os.getenv('BWA_PATH', '/usr/local/bin/bwa'),
        'samtools_path': os.getenv('SAMTOOLS_PATH', '/usr/local/bin/samtools')
    }
    
    # Initialize and run pipeline
    pipeline = VariantCallingPipeline(config)
    success = pipeline.run_pipeline(args.input, args.reference, args.output)
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main() 