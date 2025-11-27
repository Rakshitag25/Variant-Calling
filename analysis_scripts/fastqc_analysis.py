#!/usr/bin/env python3
"""
FastQC-like analysis for FASTQ files using Python
Performs comprehensive quality control analysis
"""

import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, Counter
import pyfastx
import seaborn as sns

def analyze_fastq_quality(fastq_file):
    """Perform comprehensive quality analysis on FASTQ file"""
    
    print(f"\n=== FastQC Analysis for {fastq_file} ===")
    
    # Initialize statistics
    total_reads = 0
    sequence_lengths = []
    gc_contents = []
    quality_scores = []
    per_position_quality = defaultdict(list)
    per_position_nucleotide = defaultdict(lambda: {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0})
    sequence_duplication = Counter()
    
    # Open FASTQ file
    try:
        with open(fastq_file, 'r') as f:
            lines = f.readlines()
        
        # Process in groups of 4 lines (FASTQ format)
        for i in range(0, len(lines), 4):
            if i + 3 < len(lines):
                header = lines[i].strip()
                seq = lines[i+1].strip()
                plus = lines[i+2].strip()
                qual = lines[i+3].strip()
                
                # Validate FASTQ format
                if not header.startswith('@') or plus != '+':
                    continue
                
                # Check if sequence contains only valid DNA bases
                valid_bases = set('ATCGN')
                if not all(c.upper() in valid_bases for c in seq):
                    continue
                
                # Check sequence and quality length match
                if len(seq) != len(qual):
                    continue
                
                total_reads += 1
                seq_len = len(seq)
                sequence_lengths.append(seq_len)
                
                # Calculate GC content
                gc_count = seq.upper().count('G') + seq.upper().count('C')
                gc_content = (gc_count / seq_len) * 100 if seq_len > 0 else 0
                gc_contents.append(gc_content)
                
                # Process quality scores
                qual_scores = [ord(q) - 33 for q in qual]  # Convert ASCII to Phred scores
                quality_scores.extend(qual_scores)
                
                # Per-position quality
                for pos, q_score in enumerate(qual_scores):
                    per_position_quality[pos].append(q_score)
                
                # Per-position nucleotide content
                for pos, nucleotide in enumerate(seq.upper()):
                    if nucleotide in per_position_nucleotide[pos]:
                        per_position_nucleotide[pos][nucleotide] += 1
                
                # Sequence duplication
                sequence_duplication[seq] += 1
    
    except Exception as e:
        print(f"Error reading FASTQ file: {e}")
        return None
    
    # Calculate statistics
    stats = {
        'total_sequences': total_reads,
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
    
    # Print basic statistics
    print(f"Total sequences: {stats['total_sequences']}")
    print(f"Sequence length: {stats['sequence_length']['min']}-{stats['sequence_length']['max']} bp")
    print(f"Mean sequence length: {stats['sequence_length']['mean']:.1f} bp")
    print(f"GC content: {stats['gc_content']['mean']:.1f}% ± {stats['gc_content']['std']:.1f}%")
    print(f"Quality scores: {stats['quality_scores']['min']}-{stats['quality_scores']['max']} (Phred)")
    print(f"Mean quality: {stats['quality_scores']['mean']:.1f}")
    
    # Identify potential issues
    print("\n=== Quality Assessment ===")
    
    # Check sequence length variation
    if stats['sequence_length']['min'] != stats['sequence_length']['max']:
        print("⚠️  WARNING: Variable sequence lengths detected")
    else:
        print("✅ All sequences have uniform length")
    
    # Check GC content
    if stats['gc_content']['mean'] < 20 or stats['gc_content']['mean'] > 80:
        print(f"⚠️  WARNING: Unusual GC content ({stats['gc_content']['mean']:.1f}%)")
    else:
        print(f"✅ GC content looks normal ({stats['gc_content']['mean']:.1f}%)")
    
    # Check quality scores
    if stats['quality_scores']['mean'] < 20:
        print(f"⚠️  WARNING: Low average quality score ({stats['quality_scores']['mean']:.1f})")
    elif stats['quality_scores']['mean'] >= 30:
        print(f"✅ High quality reads (Q{stats['quality_scores']['mean']:.1f})")
    else:
        print(f"✅ Good quality reads (Q{stats['quality_scores']['mean']:.1f})")
    
    # Check for duplicates
    total_unique = len(sequence_duplication)
    duplicate_rate = ((total_reads - total_unique) / total_reads) * 100 if total_reads > 0 else 0
    print(f"Sequence duplication: {duplicate_rate:.1f}% ({total_reads - total_unique}/{total_reads})")
    
    if duplicate_rate > 20:
        print("⚠️  WARNING: High duplication rate detected")
    else:
        print("✅ Duplication rate is acceptable")
    
    return stats, per_position_quality, per_position_nucleotide, gc_contents, quality_scores

def create_quality_plots(fastq_file, stats, per_position_quality, per_position_nucleotide, gc_contents, quality_scores):
    """Create FastQC-style plots"""
    
    output_dir = Path(fastq_file).parent.parent / "fastqc_output"
    output_dir.mkdir(exist_ok=True)
    
    plt.style.use('default')
    
    # 1. Per-base quality plot
    plt.figure(figsize=(12, 8))
    
    if per_position_quality:
        positions = sorted(per_position_quality.keys())
        means = [np.mean(per_position_quality[pos]) for pos in positions]
        q25 = [np.percentile(per_position_quality[pos], 25) for pos in positions]
        q75 = [np.percentile(per_position_quality[pos], 75) for pos in positions]
        
        plt.subplot(2, 2, 1)
        plt.plot(positions, means, 'b-', label='Mean Quality')
        plt.fill_between(positions, q25, q75, alpha=0.3, label='25-75 percentile')
        plt.axhline(y=30, color='g', linestyle='--', label='Q30 (99.9% accuracy)')
        plt.axhline(y=20, color='orange', linestyle='--', label='Q20 (99% accuracy)')
        plt.axhline(y=10, color='r', linestyle='--', label='Q10 (90% accuracy)')
        plt.xlabel('Position in read (bp)')
        plt.ylabel('Quality Score')
        plt.title('Per Base Quality Score')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    # 2. Per-base nucleotide content
    if per_position_nucleotide:
        plt.subplot(2, 2, 2)
        positions = sorted(per_position_nucleotide.keys())
        
        for nucleotide in ['A', 'T', 'G', 'C']:
            percentages = []
            for pos in positions:
                total = sum(per_position_nucleotide[pos].values())
                pct = (per_position_nucleotide[pos][nucleotide] / total) * 100 if total > 0 else 0
                percentages.append(pct)
            plt.plot(positions, percentages, label=f'{nucleotide}%', linewidth=2)
        
        plt.xlabel('Position in read (bp)')
        plt.ylabel('Nucleotide frequency (%)')
        plt.title('Per Base Nucleotide Content')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 100)
    
    # 3. GC content distribution
    plt.subplot(2, 2, 3)
    if gc_contents:
        plt.hist(gc_contents, bins=20, alpha=0.7, edgecolor='black')
        plt.axvline(np.mean(gc_contents), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(gc_contents):.1f}%')
        plt.xlabel('GC Content (%)')
        plt.ylabel('Number of reads')
        plt.title('GC Content Distribution')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    # 4. Quality score distribution
    plt.subplot(2, 2, 4)
    if quality_scores:
        plt.hist(quality_scores, bins=range(0, 42), alpha=0.7, edgecolor='black')
        plt.axvline(np.mean(quality_scores), color='red', linestyle='--', 
                   label=f'Mean Q: {np.mean(quality_scores):.1f}')
        plt.xlabel('Quality Score (Phred)')
        plt.ylabel('Frequency')
        plt.title('Quality Score Distribution')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    output_file = output_dir / f"{Path(fastq_file).stem}_fastqc_analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nQuality plots saved to: {output_file}")
    
    plt.show()

def main():
    # Analyze the Trial1_1.fastq file
    fastq_file = Path(r"c:\Users\PAVILION\Desktop\genome\ERR2304551_1.fastq\Trial1_1.fastq")
    
    if not fastq_file.exists():
        print(f"Error: File not found: {fastq_file}")
        return
    
    # Run analysis
    result = analyze_fastq_quality(fastq_file)
    
    if result:
        stats, per_position_quality, per_position_nucleotide, gc_contents, quality_scores = result
        
        # Create plots
        create_quality_plots(fastq_file, stats, per_position_quality, per_position_nucleotide, gc_contents, quality_scores)
        
        print(f"\n=== Analysis complete! ===")
        print("Check the generated plots for detailed quality assessment.")

if __name__ == "__main__":
    main()