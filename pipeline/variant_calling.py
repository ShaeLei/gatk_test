"""
Variant Calling Module for Variant Calling Pipeline

This module handles the core variant calling steps using GATK best practices:

1. HaplotypeCaller: Call variants using GATK's haplotype-based caller
2. VariantFiltration: Apply quality filters to variant calls
3. VariantEval: Evaluate variant call quality
4. Performance optimization for large datasets

This module demonstrates advanced variant calling techniques and optimization strategies.
"""

import logging
import os
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class HaplotypeCaller:
    """Handles variant calling using GATK HaplotypeCaller."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.gatk_path = config.get('gatk_path', 'gatk')
    
    def run_haplotype_caller(self, bam_file: str, output_dir: str) -> str:
        """
        Run GATK HaplotypeCaller for variant calling.
        
        Args:
            bam_file: Path to processed BAM file
            output_dir: Output directory for variant calls
            
        Returns:
            str: Path to raw VCF file
        """
        self.logger.info("Starting variant calling with GATK HaplotypeCaller")
        start_time = time.time()
        
        try:
            # Validate input
            if not os.path.exists(bam_file):
                raise ValueError(f"Input BAM file not found: {bam_file}")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Generate output VCF filename
            raw_vcf = os.path.join(output_dir, "raw_variants.vcf")
            
            # Optimize parameters based on BAM characteristics
            optimized_params = self._optimize_caller_parameters(bam_file)
            
            # Build GATK HaplotypeCaller command with optimized parameters
            cmd = [
                self.gatk_path, 'HaplotypeCaller',
                '--java-options', f'-Xmx{optimized_params["memory"]}g',
                '--input', bam_file,
                '--output', raw_vcf,
                '--emit-ref-confidence', 'GVCF',
                '--native-pair-hmm-threads', str(optimized_params["threads"]),
                '--max-reads-per-alignment-start', str(optimized_params["max_reads"]),
                '--min-base-quality-score', '10',
                '--min-mapping-quality', '20',
                '--standard-min-confidence-threshold-for-calling', '30.0',
                '--standard-min-confidence-threshold-for-emitting', '30.0'
            ]
            
            # Add additional parameters for performance optimization
            if optimized_params.get("use_gvcf", True):
                cmd.extend(['--emit-ref-confidence', 'GVCF'])
            
            if optimized_params.get("use_parallel", True):
                cmd.extend(['--native-pair-hmm-use-double-precision'])
            
            self.logger.info(f"Running HaplotypeCaller: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Log variant calling statistics
            self._log_variant_stats(raw_vcf)
            
            elapsed_time = time.time() - start_time
            self.logger.info(f"HaplotypeCaller completed in {elapsed_time:.2f} seconds")
            
            return raw_vcf
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"GATK HaplotypeCaller failed: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Variant calling failed: {str(e)}")
            raise
    
    def _optimize_caller_parameters(self, bam_file: str) -> Dict:
        """Optimize HaplotypeCaller parameters based on BAM characteristics."""
        self.logger.info("Optimizing HaplotypeCaller parameters")
        
        try:
            # Analyze BAM file characteristics
            file_size = os.path.getsize(bam_file) / (1024 * 1024 * 1024)  # GB
            coverage = self._estimate_coverage(bam_file)
            
            # Optimize memory allocation
            memory = self._calculate_optimal_memory(file_size, coverage)
            
            # Optimize thread allocation
            threads = self._calculate_optimal_threads()
            
            # Optimize read limits
            max_reads = self._calculate_max_reads(coverage)
            
            optimized_params = {
                'memory': memory,
                'threads': threads,
                'max_reads': max_reads,
                'use_gvcf': True,
                'use_parallel': True
            }
            
            self.logger.info(f"Optimized HaplotypeCaller parameters: {optimized_params}")
            
            return optimized_params
            
        except Exception as e:
            self.logger.error(f"Parameter optimization failed: {str(e)}")
            # Return default parameters
            return {
                'memory': self.config.get('memory', 8),
                'threads': self.config.get('threads', 4),
                'max_reads': 50,
                'use_gvcf': True,
                'use_parallel': True
            }
    
    def _estimate_coverage(self, bam_file: str) -> float:
        """Estimate average coverage from BAM file."""
        # This would use samtools to calculate coverage
        # For demonstration, returning a mock value
        return 30.0  # 30x coverage
    
    def _calculate_optimal_memory(self, file_size_gb: float, coverage: float) -> int:
        """Calculate optimal memory allocation for HaplotypeCaller."""
        # HaplotypeCaller is memory-intensive
        base_memory = max(8, int(file_size_gb * 2))  # 2GB per GB of BAM
        
        # Adjust for coverage
        if coverage > 50:
            base_memory = int(base_memory * 1.5)
        
        # Cap at available system memory
        max_memory = 24  # Leave some memory for system
        
        return min(base_memory, max_memory)
    
    def _calculate_optimal_threads(self) -> int:
        """Calculate optimal thread allocation for HaplotypeCaller."""
        # HaplotypeCaller benefits from multiple threads but not too many
        import multiprocessing
        available_cores = multiprocessing.cpu_count()
        optimal_threads = max(1, min(8, int(available_cores * 0.5)))
        
        return optimal_threads
    
    def _calculate_max_reads(self, coverage: float) -> int:
        """Calculate maximum reads per alignment start."""
        if coverage > 100:
            return 100
        elif coverage > 50:
            return 75
        else:
            return 50
    
    def _log_variant_stats(self, vcf_file: str):
        """Log variant calling statistics."""
        try:
            # Count variants in VCF file
            variant_count = self._count_variants(vcf_file)
            self.logger.info(f"HaplotypeCaller called {variant_count} variants")
            
        except Exception as e:
            self.logger.warning(f"Could not calculate variant statistics: {str(e)}")
    
    def _count_variants(self, vcf_file: str) -> int:
        """Count variants in VCF file."""
        try:
            with open(vcf_file, 'r') as f:
                count = sum(1 for line in f if not line.startswith('#'))
            return count
        except Exception:
            return 0


class VariantFiltering:
    """Handles variant filtering using GATK VariantFiltration."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.gatk_path = config.get('gatk_path', 'gatk')
    
    def run_variant_filtration(self, vcf_file: str, output_dir: str) -> str:
        """
        Run GATK VariantFiltration to apply quality filters.
        
        Args:
            vcf_file: Path to raw VCF file
            output_dir: Output directory for filtered variants
            
        Returns:
            str: Path to filtered VCF file
        """
        self.logger.info("Starting variant filtering with GATK VariantFiltration")
        start_time = time.time()
        
        try:
            # Validate input
            if not os.path.exists(vcf_file):
                raise ValueError(f"Input VCF file not found: {vcf_file}")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Generate output VCF filename
            filtered_vcf = os.path.join(output_dir, "filtered_variants.vcf")
            
            # Apply SNP filters
            snp_filtered_vcf = self._apply_snp_filters(vcf_file, output_dir)
            
            # Apply INDEL filters
            indel_filtered_vcf = self._apply_indel_filters(snp_filtered_vcf, output_dir)
            
            # Apply additional quality filters
            final_filtered_vcf = self._apply_quality_filters(indel_filtered_vcf, output_dir)
            
            # Copy final filtered VCF to output location
            import shutil
            shutil.copy2(final_filtered_vcf, filtered_vcf)
            
            # Log filtering statistics
            self._log_filtering_stats(vcf_file, filtered_vcf)
            
            elapsed_time = time.time() - start_time
            self.logger.info(f"Variant filtering completed in {elapsed_time:.2f} seconds")
            
            return filtered_vcf
            
        except Exception as e:
            self.logger.error(f"Variant filtering failed: {str(e)}")
            raise
    
    def _apply_snp_filters(self, vcf_file: str, output_dir: str) -> str:
        """Apply SNP-specific filters."""
        snp_filtered_vcf = os.path.join(output_dir, "snp_filtered.vcf")
        
        cmd = [
            self.gatk_path, 'VariantFiltration',
            '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
            '--input', vcf_file,
            '--output', snp_filtered_vcf,
            '--filter-name', 'QD2',
            '--filter-expression', 'QD < 2.0',
            '--filter-name', 'FS60',
            '--filter-expression', 'FS > 60.0',
            '--filter-name', 'SOR3',
            '--filter-expression', 'SOR > 3.0',
            '--filter-name', 'MQ40',
            '--filter-expression', 'MQ < 40.0',
            '--filter-name', 'MQRankSum-12.5',
            '--filter-expression', 'MQRankSum < -12.5',
            '--filter-name', 'ReadPosRankSum-8',
            '--filter-expression', 'ReadPosRankSum < -8.0'
        ]
        
        self.logger.info(f"Applying SNP filters: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        return snp_filtered_vcf
    
    def _apply_indel_filters(self, vcf_file: str, output_dir: str) -> str:
        """Apply INDEL-specific filters."""
        indel_filtered_vcf = os.path.join(output_dir, "indel_filtered.vcf")
        
        cmd = [
            self.gatk_path, 'VariantFiltration',
            '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
            '--input', vcf_file,
            '--output', indel_filtered_vcf,
            '--filter-name', 'QD2',
            '--filter-expression', 'QD < 2.0',
            '--filter-name', 'FS200',
            '--filter-expression', 'FS > 200.0',
            '--filter-name', 'SOR10',
            '--filter-expression', 'SOR > 10.0',
            '--filter-name', 'ReadPosRankSum-20',
            '--filter-expression', 'ReadPosRankSum < -20.0',
            '--filter-name', 'InbreedingCoeff-0.8',
            '--filter-expression', 'InbreedingCoeff < -0.8'
        ]
        
        self.logger.info(f"Applying INDEL filters: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        return indel_filtered_vcf
    
    def _apply_quality_filters(self, vcf_file: str, output_dir: str) -> str:
        """Apply additional quality filters."""
        quality_filtered_vcf = os.path.join(output_dir, "quality_filtered.vcf")
        
        cmd = [
            self.gatk_path, 'VariantFiltration',
            '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
            '--input', vcf_file,
            '--output', quality_filtered_vcf,
            '--filter-name', 'LowQual',
            '--filter-expression', 'QUAL < 30.0',
            '--filter-name', 'LowDepth',
            '--filter-expression', 'DP < 10',
            '--filter-name', 'HighDepth',
            '--filter-expression', 'DP > 1000'
        ]
        
        self.logger.info(f"Applying quality filters: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        return quality_filtered_vcf
    
    def _log_filtering_stats(self, input_vcf: str, output_vcf: str):
        """Log filtering statistics."""
        try:
            input_count = self._count_variants(input_vcf)
            output_count = self._count_variants(output_vcf)
            filtered_count = input_count - output_count
            
            self.logger.info(f"Variant filtering: {input_count} input variants, "
                           f"{output_count} output variants, {filtered_count} filtered")
            
        except Exception as e:
            self.logger.warning(f"Could not calculate filtering statistics: {str(e)}")
    
    def _count_variants(self, vcf_file: str) -> int:
        """Count variants in VCF file."""
        try:
            with open(vcf_file, 'r') as f:
                count = sum(1 for line in f if not line.startswith('#'))
            return count
        except Exception:
            return 0


class VariantEval:
    """Handles variant evaluation and quality assessment."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.gatk_path = config.get('gatk_path', 'gatk')
    
    def run_variant_eval(self, vcf_file: str, reference: str, output_dir: str) -> Dict:
        """
        Run GATK VariantEval for comprehensive variant evaluation.
        
        Args:
            vcf_file: Path to filtered VCF file
            reference: Path to reference genome
            output_dir: Output directory for evaluation results
            
        Returns:
            Dict: Evaluation metrics
        """
        self.logger.info("Starting variant evaluation with GATK VariantEval")
        
        try:
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Generate evaluation report
            eval_report = os.path.join(output_dir, "variant_eval.txt")
            
            cmd = [
                self.gatk_path, 'VariantEval',
                '--java-options', f'-Xmx{self.config.get("memory", 8)}g',
                '--input', vcf_file,
                '--reference', reference,
                '--output', eval_report,
                '--eval', 'eval:input'
            ]
            
            # Add known sites for comparison if available
            known_sites = self._get_known_sites()
            for site in known_sites:
                cmd.extend(['--comp', f'known:{site}'])
            
            self.logger.info(f"Running VariantEval: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse evaluation metrics
            metrics = self._parse_eval_metrics(eval_report)
            
            return metrics
            
        except Exception as e:
            self.logger.error(f"Variant evaluation failed: {str(e)}")
            raise
    
    def _get_known_sites(self) -> List[str]:
        """Get known variant sites for evaluation."""
        # In a real implementation, these would be paths to known variant databases
        known_sites = []
        
        # Example known sites (these would be actual file paths)
        # known_sites = [
        #     "/path/to/dbsnp_138.hg38.vcf",
        #     "/path/to/1000G_phase1.snps.high_confidence.hg38.vcf"
        # ]
        
        return known_sites
    
    def _parse_eval_metrics(self, eval_report: str) -> Dict:
        """Parse evaluation metrics from VariantEval output."""
        metrics = {
            'total_variants': 0,
            'snp_count': 0,
            'indel_count': 0,
            'transition_count': 0,
            'transversion_count': 0,
            'ti_tv_ratio': 0.0,
            'novel_variants': 0,
            'known_variants': 0
        }
        
        try:
            with open(eval_report, 'r') as f:
                lines = f.readlines()
            
            # Parse metrics from report
            for line in lines:
                if 'Total number of variants' in line:
                    metrics['total_variants'] = int(line.split()[-1])
                elif 'SNP count' in line:
                    metrics['snp_count'] = int(line.split()[-1])
                elif 'INDEL count' in line:
                    metrics['indel_count'] = int(line.split()[-1])
                elif 'Ti/Tv ratio' in line:
                    metrics['ti_tv_ratio'] = float(line.split()[-1])
            
        except Exception as e:
            self.logger.warning(f"Could not parse evaluation metrics: {str(e)}")
        
        return metrics


class VariantCallingOptimizer:
    """Handles performance optimization for variant calling."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def optimize_variant_calling(self, bam_file: str, reference: str) -> Dict:
        """
        Optimize variant calling parameters based on data characteristics.
        
        Args:
            bam_file: Path to BAM file
            reference: Path to reference genome
            
        Returns:
            Dict: Optimized parameters
        """
        self.logger.info("Optimizing variant calling parameters")
        
        try:
            # Analyze data characteristics
            coverage = self._estimate_coverage(bam_file)
            genome_size = self._estimate_genome_size(reference)
            read_length = self._estimate_read_length(bam_file)
            
            # Optimize parameters
            optimized_params = {
                'memory': self._calculate_optimal_memory(coverage, genome_size),
                'threads': self._calculate_optimal_threads(),
                'chunk_size': self._calculate_chunk_size(genome_size),
                'max_reads': self._calculate_max_reads(coverage, read_length),
                'confidence_threshold': self._calculate_confidence_threshold(coverage)
            }
            
            self.logger.info(f"Optimized variant calling parameters: {optimized_params}")
            
            return optimized_params
            
        except Exception as e:
            self.logger.error(f"Variant calling optimization failed: {str(e)}")
            return self._get_default_params()
    
    def _estimate_coverage(self, bam_file: str) -> float:
        """Estimate average coverage."""
        # This would calculate actual coverage
        return 30.0  # Mock value
    
    def _estimate_genome_size(self, reference: str) -> int:
        """Estimate genome size in base pairs."""
        # This would calculate actual genome size
        return 3000000000  # 3GB for human genome
    
    def _estimate_read_length(self, bam_file: str) -> int:
        """Estimate average read length."""
        # This would calculate actual read length
        return 150  # Mock value
    
    def _calculate_optimal_memory(self, coverage: float, genome_size: int) -> int:
        """Calculate optimal memory allocation."""
        # Higher coverage requires more memory
        base_memory = max(8, int(coverage / 10 * 4))
        return min(base_memory, 24)  # Cap at 24GB
    
    def _calculate_optimal_threads(self) -> int:
        """Calculate optimal thread allocation."""
        import multiprocessing
        available_cores = multiprocessing.cpu_count()
        return max(1, min(8, int(available_cores * 0.5)))
    
    def _calculate_chunk_size(self, genome_size: int) -> int:
        """Calculate optimal chunk size for processing."""
        # Larger genomes benefit from smaller chunks
        if genome_size > 2000000000:  # >2GB
            return 10000000  # 10MB chunks
        else:
            return 50000000  # 50MB chunks
    
    def _calculate_max_reads(self, coverage: float, read_length: int) -> int:
        """Calculate maximum reads per alignment start."""
        if coverage > 100:
            return 100
        elif coverage > 50:
            return 75
        else:
            return 50
    
    def _calculate_confidence_threshold(self, coverage: float) -> float:
        """Calculate optimal confidence threshold."""
        if coverage > 50:
            return 20.0  # Lower threshold for high coverage
        else:
            return 30.0  # Higher threshold for low coverage
    
    def _get_default_params(self) -> Dict:
        """Get default parameters."""
        return {
            'memory': 8,
            'threads': 4,
            'chunk_size': 10000000,
            'max_reads': 50,
            'confidence_threshold': 30.0
        } 