#!/usr/bin/env python3
"""
Batch Processing Script for FASTQ Chunks
Processes all chunks automatically and combines results
"""

import os
import sys
from pathlib import Path
import time
import subprocess
import json
from collections import defaultdict

def find_chunk_files(base_dir):
    """Find all chunk files in the directory"""
    chunk_files = []
    
    # Look for chunks in both ERR files
    for err_dir in ["ERR2304551_1.fastq", "ERR2304551_2.fastq"]:
        chunk_dir = Path(base_dir) / err_dir / "chunks"
        if chunk_dir.exists():
            chunks = list(chunk_dir.glob("*_chunk_*.fastq"))
            chunk_files.extend(sorted(chunks))
    
    return chunk_files

def process_single_chunk(chunk_file, output_dir):
    """Process a single chunk file using our FastQC analysis"""
    
    print(f"\nProcessing: {chunk_file.name}")
    start_time = time.time()
    
    try:
        # Import our analysis function
        sys.path.append(str(Path(__file__).parent))
        from fastqc_analysis import analyze_fastq_quality
        
        # Run analysis on chunk
        result = analyze_fastq_quality(chunk_file)
        
        if result:
            stats, per_position_quality, per_position_nucleotide, gc_contents, quality_scores = result
            
            # Save chunk results
            results_file = output_dir / f"{chunk_file.stem}_results.json"
            
            chunk_results = {
                'chunk_file': str(chunk_file),
                'processing_time': time.time() - start_time,
                'stats': stats,
                'gc_distribution': gc_contents,
                'quality_distribution': quality_scores
            }
            
            with open(results_file, 'w') as f:
                json.dump(chunk_results, f, indent=2)
            
            print(f"✅ Completed in {time.time() - start_time:.1f}s")
            return chunk_results
            
        else:
            print(f"❌ Failed to analyze {chunk_file.name}")
            return None
            
    except Exception as e:
        print(f"❌ Error processing {chunk_file.name}: {e}")
        return None

def combine_chunk_results(results_dir, output_file):
    """Combine all chunk results into summary statistics"""
    
    print(f"\n=== Combining Results ===")
    
    # Find all result files
    result_files = list(results_dir.glob("*_results.json"))
    
    if not result_files:
        print("No result files found!")
        return
    
    # Initialize combined statistics
    combined_stats = {
        'total_chunks_processed': len(result_files),
        'total_sequences': 0,
        'combined_gc_content': [],
        'combined_quality_scores': [],
        'chunk_summaries': [],
        'overall_stats': {}
    }
    
    print(f"Found {len(result_files)} result files to combine...")
    
    for result_file in result_files:
        try:
            with open(result_file, 'r') as f:
                chunk_data = json.load(f)
            
            # Add to combined statistics
            combined_stats['total_sequences'] += chunk_data['stats']['total_sequences']
            combined_stats['combined_gc_content'].extend(chunk_data['gc_distribution'])
            combined_stats['combined_quality_scores'].extend(chunk_data['quality_distribution'])
            
            # Store chunk summary
            combined_stats['chunk_summaries'].append({
                'chunk_file': chunk_data['chunk_file'],
                'sequences': chunk_data['stats']['total_sequences'],
                'mean_gc': chunk_data['stats']['gc_content']['mean'],
                'mean_quality': chunk_data['stats']['quality_scores']['mean'],
                'processing_time': chunk_data['processing_time']
            })
            
        except Exception as e:
            print(f"Warning: Could not process {result_file}: {e}")
    
    # Calculate overall statistics
    if combined_stats['combined_gc_content']:
        import numpy as np
        combined_stats['overall_stats'] = {
            'total_sequences': combined_stats['total_sequences'],
            'mean_gc_content': float(np.mean(combined_stats['combined_gc_content'])),
            'std_gc_content': float(np.std(combined_stats['combined_gc_content'])),
            'mean_quality': float(np.mean(combined_stats['combined_quality_scores'])),
            'min_quality': float(np.min(combined_stats['combined_quality_scores'])),
            'max_quality': float(np.max(combined_stats['combined_quality_scores']))
        }
    
    # Save combined results
    with open(output_file, 'w') as f:
        json.dump(combined_stats, f, indent=2)
    
    # Print summary
    print(f"\n=== COMBINED RESULTS SUMMARY ===")
    if 'overall_stats' in combined_stats:
        overall = combined_stats['overall_stats']
        print(f"Total sequences processed: {overall['total_sequences']:,}")
        print(f"Mean GC content: {overall['mean_gc_content']:.1f}%")
        print(f"Mean quality score: {overall['mean_quality']:.1f}")
        print(f"Quality range: {overall['min_quality']:.1f} - {overall['max_quality']:.1f}")
    
    print(f"Combined results saved to: {output_file}")
    
    return combined_stats

def main():
    """Main batch processing function"""
    
    base_dir = Path(r"c:\Users\PAVILION\Desktop\genome")
    output_dir = base_dir / "chunk_results"
    output_dir.mkdir(exist_ok=True)
    
    print("=== FASTQ Chunk Batch Processing ===")
    
    # Find all chunk files
    chunk_files = find_chunk_files(base_dir)
    
    if not chunk_files:
        print("No chunk files found!")
        print("Make sure you've run the chunking script first.")
        return
    
    print(f"Found {len(chunk_files)} chunk files to process")
    
    # Ask for confirmation
    print("\nChunk files to process:")
    for i, chunk_file in enumerate(chunk_files[:5]):  # Show first 5
        file_size = chunk_file.stat().st_size / (1024 * 1024)
        print(f"  {i+1}. {chunk_file.name} ({file_size:.1f} MB)")
    
    if len(chunk_files) > 5:
        print(f"  ... and {len(chunk_files) - 5} more chunks")
    
    user_input = input(f"\nProcess all {len(chunk_files)} chunks? (y/n): ")
    if user_input.lower() != 'y':
        print("Processing cancelled.")
        return
    
    # Import optimized analysis
    sys.path.append(str(Path(__file__).parent))
    from optimized_fastqc import process_chunk_batch, combine_results
    
    # Process all chunks
    print(f"\nStarting optimized batch processing...")
    start_total = time.time()
    
    results = process_chunk_batch(chunk_files, output_dir)
    
    total_time = time.time() - start_total
    
    print(f"\n=== BATCH PROCESSING COMPLETE ===")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Successful: {len(results)}/{len(chunk_files)} chunks")
    
    if results:
        # Combine all results
        combined_file = output_dir / "combined_analysis_results.json"
        combined_stats = combine_results(results, combined_file)
        
        # Create summary report
        summary_file = output_dir / "processing_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("FASTQ Chunk Processing Summary\n")
            f.write("=" * 40 + "\n")
            f.write(f"Total chunks processed: {len(results)}/{len(chunk_files)}\n")
            f.write(f"Total sequences: {combined_stats['overall_statistics']['total_sequences']:,}\n")
            f.write(f"Mean GC content: {combined_stats['overall_statistics']['mean_gc_content']:.1f}%\n")
            f.write(f"Mean quality: Q{combined_stats['overall_statistics']['mean_quality']:.1f}\n")
            f.write(f"Total processing time: {total_time/60:.1f} minutes\n")
            f.write(f"Average time per chunk: {total_time/len(chunk_files):.1f} seconds\n")
            f.write(f"Results directory: {output_dir}\n")
        
        print(f"\nDetailed summary saved to: {summary_file}")
        print(f"Individual chunk results in: {output_dir}")
    
    else:
        print("No chunks were successfully processed!")

if __name__ == "__main__":
    main()