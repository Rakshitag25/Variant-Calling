#!/usr/bin/env python3
"""
Optimized FastQC Analysis for Batch Processing
Memory-efficient analysis for large chunk processing
"""

import os
import sys
from pathlib import Path
import numpy as np
from collections import defaultdict, Counter
import json
import time

def analyze_fastq_quality(fastq_file):
    """
    Memory-efficient FastQC analysis for chunk processing
    Returns statistics without storing large datasets
    """
    
    print(f"Analyzing {Path(fastq_file).name}...", end=" ")
    start_time = time.time()
    
    # Initialize counters
    total_reads = 0
    sequence_lengths = []
    gc_sum = 0
    quality_sum = 0
    quality_count = 0
    min_quality = float('inf')
    max_quality = 0
    min_length = float('inf')
    max_length = 0
    
    # Sample data for detailed analysis (memory efficient)
    gc_sample = []
    quality_sample = []
    sample_every = 100  # Sample every 100th read for detailed stats
    
    try:
        with open(fastq_file, 'r') as f:
            line_count = 0
            
            for line in f:
                line = line.strip()
                line_count += 1
                position_in_read = (line_count - 1) % 4
                
                if position_in_read == 0:  # Header
                    if not line.startswith('@'):
                        continue
                        
                elif position_in_read == 1:  # Sequence
                    seq = line
                    seq_len = len(seq)
                    
                    # Basic validation
                    valid_bases = set('ATCGN')
                    if not all(c.upper() in valid_bases for c in seq):
                        continue
                    
                    # Update length stats
                    min_length = min(min_length, seq_len)
                    max_length = max(max_length, seq_len)
                    sequence_lengths.append(seq_len)
                    
                    # Calculate GC content
                    gc_count = seq.upper().count('G') + seq.upper().count('C')
                    gc_content = (gc_count / seq_len) * 100 if seq_len > 0 else 0
                    gc_sum += gc_content
                    
                    # Sample for detailed analysis
                    if total_reads % sample_every == 0:
                        gc_sample.append(gc_content)
                        
                elif position_in_read == 2:  # Plus line
                    if line != '+':
                        continue
                        
                elif position_in_read == 3:  # Quality
                    qual = line
                    
                    # Validate length match
                    if len(qual) != len(seq):
                        continue
                    
                    # Process quality scores
                    try:
                        qual_scores = [ord(q) - 33 for q in qual]
                        
                        # Update quality stats
                        chunk_min = min(qual_scores)
                        chunk_max = max(qual_scores)
                        chunk_sum = sum(qual_scores)
                        
                        min_quality = min(min_quality, chunk_min)
                        max_quality = max(max_quality, chunk_max)
                        quality_sum += chunk_sum
                        quality_count += len(qual_scores)
                        
                        # Sample for detailed analysis
                        if total_reads % sample_every == 0:
                            quality_sample.extend(qual_scores[:10])  # Sample first 10 positions
                            
                    except Exception:
                        continue
                    
                    # Valid read completed
                    total_reads += 1
                    
                    # Progress for large chunks
                    if total_reads % 100000 == 0:
                        elapsed = time.time() - start_time
                        print(f"\n  Progress: {total_reads:,} reads ({elapsed:.1f}s)", end=" ")
    
    except Exception as e:
        print(f"Error: {e}")
        return None
    
    # Handle edge cases
    if total_reads == 0:
        print("No valid reads found!")
        return None
    
    if min_length == float('inf'):
        min_length = 0
    if min_quality == float('inf'):
        min_quality = 0
    
    # Calculate final statistics
    stats = {
        'total_sequences': total_reads,
        'sequence_length': {
            'min': int(min_length),
            'max': int(max_length),
            'mean': float(np.mean(sequence_lengths)) if sequence_lengths else 0,
        },
        'gc_content': {
            'mean': float(gc_sum / total_reads) if total_reads > 0 else 0,
            'std': float(np.std(gc_sample)) if gc_sample else 0,
        },
        'quality_scores': {
            'mean': float(quality_sum / quality_count) if quality_count > 0 else 0,
            'min': int(min_quality),
            'max': int(max_quality),
        }
    }
    
    elapsed = time.time() - start_time
    print(f"✅ {total_reads:,} reads, Q{stats['quality_scores']['mean']:.1f}, {stats['gc_content']['mean']:.1f}% GC ({elapsed:.1f}s)")
    
    return stats, {}, {}, gc_sample, quality_sample

def process_chunk_batch(chunk_files, output_dir):
    """Process multiple chunks efficiently"""
    
    results = []
    total_start = time.time()
    
    for i, chunk_file in enumerate(chunk_files):
        print(f"\n[{i+1}/{len(chunk_files)}] ", end="")
        
        result = analyze_fastq_quality(chunk_file)
        if result:
            stats, _, _, gc_sample, qual_sample = result
            
            # Save individual chunk result
            chunk_result = {
                'chunk_file': str(chunk_file),
                'chunk_number': i + 1,
                'stats': stats,
                'gc_sample': gc_sample[:1000],  # Limit sample size
                'quality_sample': qual_sample[:1000]
            }
            
            results.append(chunk_result)
            
            # Save to individual file
            result_file = output_dir / f"{chunk_file.stem}_result.json"
            with open(result_file, 'w') as f:
                json.dump(chunk_result, f, indent=2)
    
    total_time = time.time() - total_start
    print(f"\n\nBatch complete: {len(results)}/{len(chunk_files)} chunks processed in {total_time/60:.1f} minutes")
    
    return results

def combine_results(results, output_file):
    """Combine chunk results into overall statistics"""
    
    if not results:
        print("No results to combine!")
        return
    
    print(f"\nCombining results from {len(results)} chunks...")
    
    # Aggregate statistics
    total_sequences = sum(r['stats']['total_sequences'] for r in results)
    
    # Weighted averages
    total_gc_sum = sum(r['stats']['gc_content']['mean'] * r['stats']['total_sequences'] for r in results)
    total_qual_sum = sum(r['stats']['quality_scores']['mean'] * r['stats']['total_sequences'] for r in results)
    
    overall_gc = total_gc_sum / total_sequences if total_sequences > 0 else 0
    overall_qual = total_qual_sum / total_sequences if total_sequences > 0 else 0
    
    # Combined samples
    all_gc_samples = []
    all_qual_samples = []
    
    for r in results:
        all_gc_samples.extend(r['gc_sample'])
        all_qual_samples.extend(r['quality_sample'])
    
    # Final combined statistics
    combined_stats = {
        'processing_summary': {
            'chunks_processed': len(results),
            'total_sequences': total_sequences,
            'processing_date': time.strftime('%Y-%m-%d %H:%M:%S')
        },
        'overall_statistics': {
            'total_sequences': total_sequences,
            'mean_gc_content': float(overall_gc),
            'std_gc_content': float(np.std(all_gc_samples)) if all_gc_samples else 0,
            'mean_quality': float(overall_qual),
            'min_quality': float(min(r['stats']['quality_scores']['min'] for r in results)) if results else 0,
            'max_quality': float(max(r['stats']['quality_scores']['max'] for r in results)) if results else 0,
        },
        'chunk_summaries': [
            {
                'chunk_file': r['chunk_file'],
                'sequences': r['stats']['total_sequences'],
                'gc_content': r['stats']['gc_content']['mean'],
                'quality': r['stats']['quality_scores']['mean']
            }
            for r in results
        ]
    }
    
    # Save combined results
    with open(output_file, 'w') as f:
        json.dump(combined_stats, f, indent=2)
    
    # Print summary
    print(f"\n=== FINAL COMBINED RESULTS ===")
    print(f"Total sequences: {total_sequences:,}")
    print(f"Mean GC content: {overall_gc:.1f}%")
    print(f"Mean quality score: Q{overall_qual:.1f}")
    print(f"Results saved to: {output_file}")
    
    return combined_stats

if __name__ == "__main__":
    # Test with a single chunk
    if len(sys.argv) > 1:
        chunk_file = Path(sys.argv[1])
        if chunk_file.exists():
            result = analyze_fastq_quality(chunk_file)
            if result:
                stats, _, _, _, _ = result
                print(f"\nTest analysis complete for {chunk_file.name}")
        else:
            print(f"File not found: {chunk_file}")
    else:
        print("Usage: python optimized_fastqc.py <chunk_file>")
        print("Or import this module for batch processing")