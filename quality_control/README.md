# Quality Control

This folder contains FastQC analysis scripts and quality control results:

## Files:
- `fastqc_large_files.py` - Memory-efficient FastQC analysis for large files
- `optimized_fastqc.py` - Optimized FastQC implementation for chunk processing
- `file_wise_analysis.py` - Generate file-wise quality reports and visualizations

## Directories:
- `fastqc_output/` - Raw FastQC analysis results
- `chunk_results/` - Individual chunk analysis results (JSON format)

## Usage:
These scripts and results provide comprehensive quality control analysis of FASTQ sequencing data, including per-base quality scores, GC content, sequence length distribution, and other quality metrics.