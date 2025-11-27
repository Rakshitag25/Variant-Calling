#!/usr/bin/env python3
"""
FASTQ Preprocessing Script
Performs quality control and basic preprocessing on paired-end FASTQ files
"""

import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def parse_fastq(filepath):
    """Parse FASTQ file and return read information"""
    reads = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Process in groups of 4 lines
    for i in range(0, len(lines), 4):
        if i + 3 < len(lines):
            header = lines[i].strip()
            sequence = lines[i+1].strip()
            plus = lines[i+2].strip()
            quality = lines[i+3].strip()
            
            # Validate FASTQ format
            if not header.startswith('@'):
                print(f"Warning: Invalid header at line {i+1}: {header[:50]}...")
                continue
            
            if plus != '+':
                print(f"Warning: Invalid separator at line {i+3}: {plus}")
                continue
            
            # Check if sequence contains only valid DNA bases
            valid_bases = set('ATCGN')
            if not all(c.upper() in valid_bases for c in sequence):
                print(f"Warning: Invalid sequence at line {i+2}: {sequence[:50]}...")
                continue
            
            # Check sequence and quality length match
            if len(sequence) != len(quality):
                print(f"Warning: Sequence/quality length mismatch at line {i+2}: seq={len(sequence)}, qual={len(quality)}")
                continue
            
            reads.append({
                'header': header,
                'sequence': sequence,
                'quality': quality,
                'length': len(sequence)
            })
    
    return reads

def quality_stats(reads):
    """Calculate quality statistics"""
    if not reads:
        return {}
    
    lengths = [read['length'] for read in reads]
    gc_counts = []
    
    for read in reads:
        seq = read['sequence'].upper()
        gc_count = seq.count('G') + seq.count('C')
        gc_content = (gc_count / len(seq)) * 100 if len(seq) > 0 else 0
        gc_counts.append(gc_content)
    
    stats = {
        'total_reads': len(reads),
        'min_length': min(lengths) if lengths else 0,
        'max_length': max(lengths) if lengths else 0,
        'avg_length': np.mean(lengths) if lengths else 0,
        'avg_gc_content': np.mean(gc_counts) if gc_counts else 0,
        'length_distribution': lengths,
        'gc_distribution': gc_counts
    }
    
    return stats

def filter_reads(reads, min_length=50):
    """Filter reads based on quality criteria"""
    filtered = []
    
    for read in reads:
        # Basic filters
        if len(read['sequence']) >= min_length:
            # Check for excessive Ns
            n_count = read['sequence'].upper().count('N')
            n_percentage = (n_count / len(read['sequence'])) * 100
            
            if n_percentage < 10:  # Less than 10% Ns
                filtered.append(read)
    
    return filtered

def preprocess_fastq(file_path, output_dir):
    """Main preprocessing function"""
    print(f"\n=== Processing {file_path} ===")
    
    # Parse FASTQ
    reads = parse_fastq(file_path)
    print(f"Total reads found: {len(reads)}")
    
    # Calculate statistics
    stats = quality_stats(reads)
    
    print(f"Read length range: {stats['min_length']} - {stats['max_length']} bp")
    print(f"Average read length: {stats['avg_length']:.1f} bp")
    print(f"Average GC content: {stats['avg_gc_content']:.1f}%")
    
    # Filter reads
    filtered_reads = filter_reads(reads, min_length=100)
    print(f"Reads after filtering: {len(filtered_reads)}")
    print(f"Filtered out: {len(reads) - len(filtered_reads)} reads")
    
    # Save filtered reads
    filename = Path(file_path).stem
    output_file = Path(output_dir) / f"{filename}_filtered.fastq"
    
    with open(output_file, 'w') as f:
        for read in filtered_reads:
            f.write(f"{read['header']}\n")
            f.write(f"{read['sequence']}\n")
            f.write("+\n")
            f.write(f"{read['quality']}\n")
    
    print(f"Filtered file saved: {output_file}")
    
    return stats, filtered_reads

def main():
    # Set up paths
    genome_dir = Path(r"c:\Users\PAVILION\Desktop\genome")
    output_dir = genome_dir / "preprocessed"
    output_dir.mkdir(exist_ok=True)
    
    # Find FASTQ files
    fastq_files = [
        genome_dir / "ERR2304551_1.fastq" / "Trial1_1.fastq",
    ]
    
    # Process each file
    all_stats = {}
    
    for file_path in fastq_files:
        if file_path.exists():
            try:
                stats, filtered_reads = preprocess_fastq(file_path, output_dir)
                all_stats[file_path.name] = stats
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
        else:
            print(f"File not found: {file_path}")
    
    # Print summary
    print("\n=== PREPROCESSING SUMMARY ===")
    for filename, stats in all_stats.items():
        print(f"\n{filename}:")
        print(f"  Total reads: {stats['total_reads']}")
        print(f"  Average length: {stats['avg_length']:.1f} bp")
        print(f"  GC content: {stats['avg_gc_content']:.1f}%")

if __name__ == "__main__":
    main()