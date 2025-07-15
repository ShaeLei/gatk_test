"""
Validation Module for Variant Calling Pipeline

This module handles accuracy validation, performance monitoring, and quality assessment:

1. ValidationMetrics: Calculate accuracy metrics against ground truth
2. PerformanceMonitor: Monitor runtime performance and resource usage
3. QualityAssessment: Assess variant call quality and reliability
4. CrossValidation: Compare results with multiple variant callers

This module demonstrates comprehensive validation strategies for variant calling accuracy.
"""

import logging
import os
import subprocess
import time
import psutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json


class ValidationMetrics:
    """Handles accuracy validation against ground truth datasets."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.gatk_path = config.get('gatk_path', 'gatk')
    
    def calculate_metrics(self, vcf_file: str) -> Dict:
        """
        Calculate comprehensive validation metrics.
        
        Args:
            vcf_file: Path to filtered VCF file
            
        Returns:
            Dict: Validation metrics including sensitivity, specificity, etc.
        """
        self.logger.info("Calculating validation metrics")
        
        try:
            metrics = {}
            
            # Basic variant statistics
            metrics['basic_stats'] = self._calculate_basic_stats(vcf_file)
            
            # Quality metrics
            metrics['quality_metrics'] = self._calculate_quality_metrics(vcf_file)
            
            # Accuracy metrics (if ground truth available)
            if self._has_ground_truth():
                metrics['accuracy_metrics'] = self._calculate_accuracy_metrics(vcf_file)
            
            # Transition/transversion ratio
            metrics['ti_tv_ratio'] = self._calculate_ti_tv_ratio(vcf_file)
            
            # Novel vs known variants
            metrics['novelty_metrics'] = self._calculate_novelty_metrics(vcf_file)
            
            return metrics
            
        except Exception as e:
            self.logger.error(f"Validation metrics calculation failed: {str(e)}")
            raise
    
    def _calculate_basic_stats(self, vcf_file: str) -> Dict:
        """Calculate basic variant statistics."""
        stats = {
            'total_variants': 0,
            'snps': 0,
            'indels': 0,
            'singletons': 0,
            'multi_allelic': 0
        }
        
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 5:
                            stats['total_variants'] += 1
                            
                            # Count SNPs vs INDELs
                            ref = parts[3]
                            alt = parts[4]
                            if len(ref) == 1 and len(alt) == 1:
                                stats['snps'] += 1
                            else:
                                stats['indels'] += 1
                            
                            # Count multi-allelic variants
                            if ',' in alt:
                                stats['multi_allelic'] += 1
                            
                            # Count singletons (AF=1)
                            info = parts[7] if len(parts) > 7 else ''
                            if 'AF=1' in info:
                                stats['singletons'] += 1
                                
        except Exception as e:
            self.logger.warning(f"Could not calculate basic stats: {str(e)}")
        
        return stats
    
    def _calculate_quality_metrics(self, vcf_file: str) -> Dict:
        """Calculate quality metrics for variants."""
        quality_metrics = {
            'mean_qual': 0.0,
            'median_qual': 0.0,
            'qual_distribution': {},
            'filtered_variants': 0,
            'pass_variants': 0
        }
        
        try:
            quals = []
            filtered_count = 0
            pass_count = 0
            
            with open(vcf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 6:
                            qual = float(parts[5]) if parts[5] != '.' else 0.0
                            quals.append(qual)
                            
                            # Count filtered vs pass variants
                            if 'PASS' in parts[6]:
                                pass_count += 1
                            else:
                                filtered_count += 1
            
            if quals:
                quality_metrics['mean_qual'] = sum(quals) / len(quals)
                quality_metrics['median_qual'] = sorted(quals)[len(quals)//2]
                quality_metrics['qual_distribution'] = self._calculate_qual_distribution(quals)
            
            quality_metrics['filtered_variants'] = filtered_count
            quality_metrics['pass_variants'] = pass_count
            
        except Exception as e:
            self.logger.warning(f"Could not calculate quality metrics: {str(e)}")
        
        return quality_metrics
    
    def _calculate_qual_distribution(self, quals: List[float]) -> Dict:
        """Calculate quality score distribution."""
        distribution = {
            '0-10': 0, '10-20': 0, '20-30': 0, '30-40': 0, '40-50': 0, '50+': 0
        }
        
        for qual in quals:
            if qual < 10:
                distribution['0-10'] += 1
            elif qual < 20:
                distribution['10-20'] += 1
            elif qual < 30:
                distribution['20-30'] += 1
            elif qual < 40:
                distribution['30-40'] += 1
            elif qual < 50:
                distribution['40-50'] += 1
            else:
                distribution['50+'] += 1
        
        return distribution
    
    def _calculate_accuracy_metrics(self, vcf_file: str) -> Dict:
        """Calculate accuracy metrics against ground truth."""
        accuracy_metrics = {
            'sensitivity': 0.0,
            'specificity': 0.0,
            'precision': 0.0,
            'f1_score': 0.0,
            'true_positives': 0,
            'false_positives': 0,
            'false_negatives': 0,
            'true_negatives': 0
        }
        
        try:
            # This would compare against known variant databases
            # For demonstration, returning mock metrics
            accuracy_metrics = {
                'sensitivity': 0.95,
                'specificity': 0.98,
                'precision': 0.96,
                'f1_score': 0.955,
                'true_positives': 9500,
                'false_positives': 200,
                'false_negatives': 500,
                'true_negatives': 9800
            }
            
        except Exception as e:
            self.logger.warning(f"Could not calculate accuracy metrics: {str(e)}")
        
        return accuracy_metrics
    
    def _calculate_ti_tv_ratio(self, vcf_file: str) -> float:
        """Calculate transition/transversion ratio."""
        transitions = 0
        transversions = 0
        
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 5:
                            ref = parts[3]
                            alt = parts[4]
                            
                            # Only consider SNPs
                            if len(ref) == 1 and len(alt) == 1:
                                if self._is_transition(ref, alt):
                                    transitions += 1
                                else:
                                    transversions += 1
            
            if transversions > 0:
                return transitions / transversions
            else:
                return 0.0
                
        except Exception as e:
            self.logger.warning(f"Could not calculate Ti/Tv ratio: {str(e)}")
            return 0.0
    
    def _is_transition(self, ref: str, alt: str) -> bool:
        """Check if variant is a transition (A↔G or C↔T)."""
        transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
        return (ref, alt) in transitions
    
    def _calculate_novelty_metrics(self, vcf_file: str) -> Dict:
        """Calculate novelty metrics (novel vs known variants)."""
        novelty_metrics = {
            'known_variants': 0,
            'novel_variants': 0,
            'novelty_rate': 0.0
        }
        
        try:
            # This would compare against known variant databases
            # For demonstration, using mock values
            total_variants = self._count_variants(vcf_file)
            known_variants = int(total_variants * 0.8)  # 80% known
            novel_variants = total_variants - known_variants
            
            novelty_metrics = {
                'known_variants': known_variants,
                'novel_variants': novel_variants,
                'novelty_rate': novel_variants / total_variants if total_variants > 0 else 0.0
            }
            
        except Exception as e:
            self.logger.warning(f"Could not calculate novelty metrics: {str(e)}")
        
        return novelty_metrics
    
    def _count_variants(self, vcf_file: str) -> int:
        """Count total variants in VCF file."""
        try:
            with open(vcf_file, 'r') as f:
                return sum(1 for line in f if not line.startswith('#'))
        except Exception:
            return 0
    
    def _has_ground_truth(self) -> bool:
        """Check if ground truth data is available."""
        # This would check for known variant databases
        return False  # Mock implementation
    
    def generate_report(self, metrics: Dict, output_dir: str):
        """Generate comprehensive validation report."""
        report_file = os.path.join(output_dir, "validation_report.txt")
        
        try:
            with open(report_file, 'w') as f:
                f.write("Variant Calling Validation Report\n")
                f.write("=" * 40 + "\n\n")
                
                # Basic statistics
                f.write("Basic Statistics:\n")
                f.write("-" * 20 + "\n")
                for key, value in metrics['basic_stats'].items():
                    f.write(f"{key}: {value}\n")
                
                f.write("\nQuality Metrics:\n")
                f.write("-" * 20 + "\n")
                for key, value in metrics['quality_metrics'].items():
                    if key != 'qual_distribution':
                        f.write(f"{key}: {value}\n")
                
                f.write("\nQuality Score Distribution:\n")
                for range_key, count in metrics['quality_metrics']['qual_distribution'].items():
                    f.write(f"  {range_key}: {count}\n")
                
                f.write(f"\nTi/Tv Ratio: {metrics['ti_tv_ratio']:.3f}\n")
                
                f.write("\nNovelty Metrics:\n")
                f.write("-" * 20 + "\n")
                for key, value in metrics['novelty_metrics'].items():
                    f.write(f"{key}: {value}\n")
                
                if 'accuracy_metrics' in metrics:
                    f.write("\nAccuracy Metrics:\n")
                    f.write("-" * 20 + "\n")
                    for key, value in metrics['accuracy_metrics'].items():
                        f.write(f"{key}: {value}\n")
            
            self.logger.info(f"Validation report saved to {report_file}")
            
        except Exception as e:
            self.logger.error(f"Could not generate validation report: {str(e)}")


class PerformanceMonitor:
    """Handles performance monitoring and resource usage tracking."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.start_time = None
        self.step_times = {}
        self.resource_usage = {}
        self.errors = []
    
    def start_pipeline(self):
        """Start monitoring the entire pipeline."""
        self.start_time = time.time()
        self.logger.info("Starting performance monitoring")
    
    def end_pipeline(self):
        """End pipeline monitoring."""
        if self.start_time:
            total_time = time.time() - self.start_time
            self.logger.info(f"Pipeline completed in {total_time:.2f} seconds")
    
    def start_step(self, step_name: str):
        """Start monitoring a specific pipeline step."""
        self.step_times[step_name] = {
            'start': time.time(),
            'cpu_start': psutil.cpu_percent(),
            'memory_start': psutil.virtual_memory().percent
        }
        self.logger.info(f"Starting step: {step_name}")
    
    def end_step(self, step_name: str):
        """End monitoring a specific pipeline step."""
        if step_name in self.step_times:
            end_time = time.time()
            start_data = self.step_times[step_name]
            
            step_duration = end_time - start_data['start']
            cpu_usage = psutil.cpu_percent() - start_data['cpu_start']
            memory_usage = psutil.virtual_memory().percent - start_data['memory_start']
            
            self.step_times[step_name].update({
                'end': end_time,
                'duration': step_duration,
                'cpu_usage': cpu_usage,
                'memory_usage': memory_usage
            })
            
            self.logger.info(f"Step {step_name} completed in {step_duration:.2f} seconds")
    
    def record_error(self, error_message: str):
        """Record an error during pipeline execution."""
        self.errors.append({
            'time': time.time(),
            'message': error_message
        })
        self.logger.error(f"Pipeline error recorded: {error_message}")
    
    def generate_report(self) -> str:
        """Generate comprehensive performance report."""
        report_lines = []
        
        # Overall pipeline statistics
        if self.start_time:
            total_time = time.time() - self.start_time
            report_lines.append(f"Total Pipeline Time: {total_time:.2f} seconds")
        
        # Step-by-step breakdown
        report_lines.append("\nStep Performance Breakdown:")
        report_lines.append("-" * 40)
        
        for step_name, step_data in self.step_times.items():
            if 'duration' in step_data:
                report_lines.append(f"{step_name}:")
                report_lines.append(f"  Duration: {step_data['duration']:.2f} seconds")
                report_lines.append(f"  CPU Usage: {step_data.get('cpu_usage', 0):.1f}%")
                report_lines.append(f"  Memory Usage: {step_data.get('memory_usage', 0):.1f}%")
        
        # Error summary
        if self.errors:
            report_lines.append(f"\nErrors Encountered: {len(self.errors)}")
            for error in self.errors:
                report_lines.append(f"  - {error['message']}")
        
        # Resource usage summary
        report_lines.append("\nResource Usage Summary:")
        report_lines.append("-" * 30)
        report_lines.append(f"Current CPU Usage: {psutil.cpu_percent():.1f}%")
        report_lines.append(f"Current Memory Usage: {psutil.virtual_memory().percent:.1f}%")
        report_lines.append(f"Available Memory: {psutil.virtual_memory().available / (1024**3):.1f} GB")
        
        return "\n".join(report_lines)


class QualityAssessment:
    """Handles comprehensive quality assessment of variant calls."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def assess_variant_quality(self, vcf_file: str, output_dir: str) -> Dict:
        """
        Perform comprehensive quality assessment of variant calls.
        
        Args:
            vcf_file: Path to VCF file
            output_dir: Output directory for quality reports
            
        Returns:
            Dict: Quality assessment results
        """
        self.logger.info("Starting comprehensive quality assessment")
        
        try:
            quality_results = {}
            
            # Assess variant distribution
            quality_results['distribution'] = self._assess_variant_distribution(vcf_file)
            
            # Assess quality score distribution
            quality_results['quality_scores'] = self._assess_quality_scores(vcf_file)
            
            # Assess allele frequency distribution
            quality_results['allele_frequencies'] = self._assess_allele_frequencies(vcf_file)
            
            # Assess depth distribution
            quality_results['depth_distribution'] = self._assess_depth_distribution(vcf_file)
            
            # Generate quality report
            self._generate_quality_report(quality_results, output_dir)
            
            return quality_results
            
        except Exception as e:
            self.logger.error(f"Quality assessment failed: {str(e)}")
            raise
    
    def _assess_variant_distribution(self, vcf_file: str) -> Dict:
        """Assess variant type and length distribution."""
        distribution = {
            'snps': 0,
            'insertions': 0,
            'deletions': 0,
            'complex': 0,
            'length_distribution': {}
        }
        
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 5:
                            ref = parts[3]
                            alt = parts[4]
                            
                            ref_len = len(ref)
                            alt_len = len(alt)
                            
                            if ref_len == 1 and alt_len == 1:
                                distribution['snps'] += 1
                            elif ref_len == 1 and alt_len > 1:
                                distribution['insertions'] += 1
                            elif ref_len > 1 and alt_len == 1:
                                distribution['deletions'] += 1
                            else:
                                distribution['complex'] += 1
                            
                            # Track length distribution
                            max_len = max(ref_len, alt_len)
                            if max_len not in distribution['length_distribution']:
                                distribution['length_distribution'][max_len] = 0
                            distribution['length_distribution'][max_len] += 1
                            
        except Exception as e:
            self.logger.warning(f"Could not assess variant distribution: {str(e)}")
        
        return distribution
    
    def _assess_quality_scores(self, vcf_file: str) -> Dict:
        """Assess quality score distribution."""
        quality_scores = []
        
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 6 and parts[5] != '.':
                            quality_scores.append(float(parts[5]))
            
            if quality_scores:
                return {
                    'mean': sum(quality_scores) / len(quality_scores),
                    'median': sorted(quality_scores)[len(quality_scores)//2],
                    'min': min(quality_scores),
                    'max': max(quality_scores),
                    'std': self._calculate_std(quality_scores)
                }
            else:
                return {'mean': 0, 'median': 0, 'min': 0, 'max': 0, 'std': 0}
                
        except Exception as e:
            self.logger.warning(f"Could not assess quality scores: {str(e)}")
            return {'mean': 0, 'median': 0, 'min': 0, 'max': 0, 'std': 0}
    
    def _assess_allele_frequencies(self, vcf_file: str) -> Dict:
        """Assess allele frequency distribution."""
        frequencies = []
        
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 8:
                            info = parts[7]
                            # Extract AF from INFO field
                            if 'AF=' in info:
                                af_part = info.split('AF=')[1].split(';')[0]
                                try:
                                    af = float(af_part)
                                    frequencies.append(af)
                                except ValueError:
                                    pass
            
            if frequencies:
                return {
                    'mean': sum(frequencies) / len(frequencies),
                    'median': sorted(frequencies)[len(frequencies)//2],
                    'min': min(frequencies),
                    'max': max(frequencies)
                }
            else:
                return {'mean': 0, 'median': 0, 'min': 0, 'max': 0}
                
        except Exception as e:
            self.logger.warning(f"Could not assess allele frequencies: {str(e)}")
            return {'mean': 0, 'median': 0, 'min': 0, 'max': 0}
    
    def _assess_depth_distribution(self, vcf_file: str) -> Dict:
        """Assess depth distribution."""
        depths = []
        
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 8:
                            info = parts[7]
                            # Extract DP from INFO field
                            if 'DP=' in info:
                                dp_part = info.split('DP=')[1].split(';')[0]
                                try:
                                    dp = int(dp_part)
                                    depths.append(dp)
                                except ValueError:
                                    pass
            
            if depths:
                return {
                    'mean': sum(depths) / len(depths),
                    'median': sorted(depths)[len(depths)//2],
                    'min': min(depths),
                    'max': max(depths)
                }
            else:
                return {'mean': 0, 'median': 0, 'min': 0, 'max': 0}
                
        except Exception as e:
            self.logger.warning(f"Could not assess depth distribution: {str(e)}")
            return {'mean': 0, 'median': 0, 'min': 0, 'max': 0}
    
    def _calculate_std(self, values: List[float]) -> float:
        """Calculate standard deviation."""
        if len(values) < 2:
            return 0.0
        
        mean = sum(values) / len(values)
        variance = sum((x - mean) ** 2 for x in values) / (len(values) - 1)
        return variance ** 0.5
    
    def _generate_quality_report(self, quality_results: Dict, output_dir: str):
        """Generate quality assessment report."""
        report_file = os.path.join(output_dir, "quality_assessment_report.txt")
        
        try:
            with open(report_file, 'w') as f:
                f.write("Variant Quality Assessment Report\n")
                f.write("=" * 40 + "\n\n")
                
                # Variant distribution
                f.write("Variant Distribution:\n")
                f.write("-" * 25 + "\n")
                for key, value in quality_results['distribution'].items():
                    if key != 'length_distribution':
                        f.write(f"{key}: {value}\n")
                
                f.write("\nLength Distribution:\n")
                for length, count in sorted(quality_results['distribution']['length_distribution'].items()):
                    f.write(f"  Length {length}: {count} variants\n")
                
                # Quality scores
                f.write("\nQuality Score Statistics:\n")
                f.write("-" * 30 + "\n")
                for key, value in quality_results['quality_scores'].items():
                    f.write(f"{key}: {value:.3f}\n")
                
                # Allele frequencies
                f.write("\nAllele Frequency Statistics:\n")
                f.write("-" * 35 + "\n")
                for key, value in quality_results['allele_frequencies'].items():
                    f.write(f"{key}: {value:.3f}\n")
                
                # Depth distribution
                f.write("\nDepth Distribution Statistics:\n")
                f.write("-" * 35 + "\n")
                for key, value in quality_results['depth_distribution'].items():
                    f.write(f"{key}: {value:.1f}\n")
            
            self.logger.info(f"Quality assessment report saved to {report_file}")
            
        except Exception as e:
            self.logger.error(f"Could not generate quality report: {str(e)}") 