"""
Variant Calling Pipeline Package

A comprehensive variant calling pipeline following GATK best practices,
optimized for accuracy and runtime performance.

This package provides:
- Complete preprocessing workflow
- BWA alignment and SAM processing
- GATK post-processing steps
- Variant calling and filtering
- Comprehensive validation and quality assessment
- Performance monitoring and optimization
"""

__version__ = "1.0.0"
__author__ = "Bioinformatics Team"
__description__ = "GATK Best Practices Variant Calling Pipeline"

# Import main pipeline components
from .main import VariantCallingPipeline
from .preprocessing import QualityControl, AdapterTrimming, QualityTrimming
from .alignment import BWAAlignment, SAMProcessing
from .postprocessing import Deduplication, BaseRecalibration
from .variant_calling import HaplotypeCaller, VariantFiltering
from .validation import ValidationMetrics, PerformanceMonitor

__all__ = [
    'VariantCallingPipeline',
    'QualityControl',
    'AdapterTrimming', 
    'QualityTrimming',
    'BWAAlignment',
    'SAMProcessing',
    'Deduplication',
    'BaseRecalibration',
    'HaplotypeCaller',
    'VariantFiltering',
    'ValidationMetrics',
    'PerformanceMonitor'
] 