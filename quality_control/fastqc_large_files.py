#!/usr/bin/env python3
"""
FastQC Analysis for Large FASTQ Files
Processes files in streaming mode to handle large datasets efficiently
"""

import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, Counter
import gc  # Garbage collection

def analyze_large_fastq_streaming(fastq_file, sample_size=10000):
    """
    Stream-based analysis for large FASTQ files
    Processes only a sample of reads to generate statistics efficiently
    """
    
    print(f"\n=== Streaming FastQC Analysis for {fastq_file} ===")
    print(f"Sampling every {sample_size//1000}k reads for analysis...")
    
    # Initialize statistics
    total_reads = 0
    analyzed_reads = 0
    sequence_lengths = []
    gc_contents = []
    quality_scores = []
    per_position_quality = defaultdict(list)
    per_position_nucleotide = defaultdict(lambda: {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0})
    
    try:
        with open(fastq_file, 'r') as f:
            line_count = 0
            current_read = {}
            
            for line in f:
                line = line.strip()
                line_count += 1
                
                # Parse FASTQ in groups of 4 lines
                position_in_read = (line_count - 1) % 4
                
                if position_in_read == 0:  # Header line
                    if line.startswith('@'):
                        current_read['header'] = line
                        total_reads += 1
                        
                        # Sample every nth read for detailed analysis
                        if total_reads % (sample_size // 1000) == 0:
                            current_read['analyze'] = True
                            analyzed_reads += 1
                        else:
                            current_read['analyze'] = False
                            
                elif position_in_read == 1:  # Sequence line
                    current_read['sequence'] = line
                    
                elif position_in_read == 2:  # Plus line
                    current_read['plus'] = line
                    
                elif position_in_read == 3:  # Quality line
                    current_read['quality'] = line
                    
                    # Process the complete read
                    if current_read.get('analyze', False):
                        process_read_for_analysis(current_read, sequence_lengths, gc_contents, 
                                                quality_scores, per_position_quality, 
                                                per_position_nucleotide)
                    
                    # Clear memory
                    current_read = {}
                    
                    # Progress indicator
                    if total_reads % 50000 == 0:
                        print(f"Processed {total_reads:,} reads...")
                        
                    # Optional: Limit total processing for very large files
                    if total_reads >= sample_size:
                        print(f"Reached sample limit of {sample_size:,} reads")
                        break
    
    except Exception as e:
        print(f"Error processing file: {e}")
        return None
    
    print(f"Analysis complete: {total_reads:,} total reads, {analyzed_reads:,} analyzed in detail")
    
    # Calculate final statistics
    stats = calculate_statistics(sequence_lengths, gc_contents, quality_scores, total_reads, analyzed_reads)
    
    return stats, per_position_quality, per_position_nucleotide, gc_contents, quality_scores

def process_read_for_analysis(read, sequence_lengths, gc_contents, quality_scores, 
                            per_position_quality, per_position_nucleotide):
    """Process a single read for detailed analysis"""
    
    seq = read['sequence']
    qual = read['quality']
    
    # Basic validation
    if len(seq) != len(qual):
        return
    
    valid_bases = set('ATCGN')
    if not all(c.upper() in valid_bases for c in seq):
        return
    
    # Sequence length
    seq_len = len(seq)
    sequence_lengths.append(seq_len)
    
    # GC content
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    gc_content = (gc_count / seq_len) * 100 if seq_len > 0 else 0
    gc_contents.append(gc_content)
    
    # Quality scores
    try:
        qual_scores = [ord(q) - 33 for q in qual]  # Convert ASCII to Phred scores
        quality_scores.extend(qual_scores)
        
        # Per-position quality (sample first 100 positions to save memory)
        for pos, q_score in enumerate(qual_scores[:100]):
            per_position_quality[pos].append(q_score)
        
        # Per-position nucleotide content (sample first 100 positions)
        for pos, nucleotide in enumerate(seq.upper()[:100]):
            if nucleotide in per_position_nucleotide[pos]:
                per_position_nucleotide[pos][nucleotide] += 1
                
    except Exception as e:
        print(f"Warning: Quality score processing error: {e}")

def calculate_statistics(sequence_lengths, gc_contents, quality_scores, total_reads, analyzed_reads):
    """Calculate final statistics"""
    
    stats = {
        'total_sequences': total_reads,
        'analyzed_sequences': analyzed_reads,
        'sequence_length': {
            'min': min(sequence_lengths) if sequence_lengths else 0,
            'max': max(sequence_lengths) if sequence_lengths else 0,
            'mean': np.mean(sequence_lengths) if sequence_lengths else 0,
        },
        'gc_content': {
            'mean': np.mean(gc_contents) if gc_contents else 0,
            'std': np.std(gc_contents) if gc_contents else 0,
        },
        'quality_scores': {
            'mean': np.mean(quality_scores) if quality_scores else 0,
            'min': min(quality_scores) if quality_scores else 0,
            'max': max(quality_scores) if quality_scores else 0,
        }
    }
    
    # Print statistics
    print(f"\n=== Results ===")
    print(f"Total sequences: {stats['total_sequences']:,}")
    print(f"Analyzed sequences: {stats['analyzed_sequences']:,}")
    print(f"Sequence length: {stats['sequence_length']['min']}-{stats['sequence_length']['max']} bp")
    print(f"Mean sequence length: {stats['sequence_length']['mean']:.1f} bp")
    print(f"GC content: {stats['gc_content']['mean']:.1f}% ± {stats['gc_content']['std']:.1f}%")
    print(f"Quality scores: {stats['quality_scores']['min']}-{stats['quality_scores']['max']} (Phred)")
    print(f"Mean quality: {stats['quality_scores']['mean']:.1f}")
    
    return stats

def analyze_file_size_and_estimate(fastq_file):
    """Estimate file processing requirements"""
    
    file_size = os.path.getsize(fastq_file)
    file_size_mb = file_size / (1024 * 1024)
    file_size_gb = file_size_mb / 1024
    
    print(f"\n=== File Analysis ===")
    print(f"File: {fastq_file}")
    print(f"File size: {file_size_mb:.1f} MB ({file_size_gb:.2f} GB)")
    
    # Estimate number of reads (assuming ~100 bytes per read on average)
    estimated_reads = file_size // 100
    print(f"Estimated reads: {estimated_reads:,}")
    
    # Memory estimation
    if file_size_gb > 2:
        print("⚠️  Large file detected - will use streaming analysis")
        recommended_sample = min(100000, estimated_reads // 10)
    else:
        print("✅ Medium file - can process with sampling")
        recommended_sample = min(50000, estimated_reads // 5)
    
    print(f"Recommended sample size: {recommended_sample:,}")
    
    return recommended_sample

def main():
    """Main function to process large FASTQ files"""
    
    # List of files to process
    fastq_files = [
        Path(r"c:\Users\PAVILION\Desktop\genome\ERR2304551_1.fastq\ERR2304551_1.fastq"),
        Path(r"c:\Users\PAVILION\Desktop\genome\ERR2304551_2.fastq\ERR2304551_2.fastq"),
    ]
    
    output_dir = Path(r"c:\Users\PAVILION\Desktop\genome\fastqc_output")
    output_dir.mkdir(exist_ok=True)
    
    for fastq_file in fastq_files:
        if fastq_file.exists():
            print(f"\n{'='*60}")
            
            # Analyze file size and get recommendations
            recommended_sample = analyze_file_size_and_estimate(fastq_file)
            
            # Run streaming analysis
            result = analyze_large_fastq_streaming(fastq_file, sample_size=recommended_sample)
            
            if result:
                stats, per_position_quality, per_position_nucleotide, gc_contents, quality_scores = result
                
                # Save results to file
                results_file = output_dir / f"{fastq_file.stem}_fastqc_results.txt"
                with open(results_file, 'w') as f:
                    f.write(f"FastQC Analysis Results for {fastq_file.name}\n")
                    f.write("="*50 + "\n")
                    f.write(f"Total sequences: {stats['total_sequences']:,}\n")
                    f.write(f"Analyzed sequences: {stats['analyzed_sequences']:,}\n")
                    f.write(f"Mean length: {stats['sequence_length']['mean']:.1f} bp\n")
                    f.write(f"GC content: {stats['gc_content']['mean']:.1f}%\n")
                    f.write(f"Mean quality: {stats['quality_scores']['mean']:.1f}\n")
                
                print(f"Results saved to: {results_file}")
                
                # Force garbage collection
                gc.collect()
            
        else:
            print(f"File not found: {fastq_file}")

if __name__ == "__main__":
    main()