#!/usr/bin/env python3
"""
File-wise Results Analysis and Visualization
Combines chunk results by original file and creates comprehensive plots
"""

import os
import sys
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import defaultdict

def load_chunk_results(results_dir):
    """Load all chunk results and organize by file"""
    
    print("Loading chunk results...")
    results_dir = Path(results_dir)
    
    # Organize results by original file
    file_results = defaultdict(list)
    
    # Find all chunk result files
    chunk_files = list(results_dir.glob("*_chunk_*_result.json"))
    
    for chunk_file in chunk_files:
        try:
            with open(chunk_file, 'r') as f:
                chunk_data = json.load(f)
            
            # Extract original file name from chunk filename
            chunk_name = chunk_file.stem.replace('_result', '')
            
            if 'ERR2304551_1' in chunk_name:
                original_file = 'ERR2304551_1.fastq'
            elif 'ERR2304551_2' in chunk_name:
                original_file = 'ERR2304551_2.fastq'
            else:
                continue
            
            file_results[original_file].append(chunk_data)
            
        except Exception as e:
            print(f"Warning: Could not load {chunk_file}: {e}")
    
    print(f"Loaded results for {len(file_results)} files")
    for file_name, chunks in file_results.items():
        print(f"  {file_name}: {len(chunks)} chunks")
    
    return file_results

def combine_file_statistics(chunks):
    """Combine statistics from all chunks of a single file"""
    
    if not chunks:
        return None
    
    # Aggregate basic stats
    total_sequences = sum(chunk['stats']['total_sequences'] for chunk in chunks)
    
    # Weighted averages for GC content and quality
    total_gc_sum = sum(chunk['stats']['gc_content']['mean'] * chunk['stats']['total_sequences'] for chunk in chunks)
    total_qual_sum = sum(chunk['stats']['quality_scores']['mean'] * chunk['stats']['total_sequences'] for chunk in chunks)
    
    mean_gc = total_gc_sum / total_sequences if total_sequences > 0 else 0
    mean_quality = total_qual_sum / total_sequences if total_sequences > 0 else 0
    
    # Collect all sample data
    all_gc_samples = []
    all_quality_samples = []
    
    for chunk in chunks:
        if 'gc_sample' in chunk:
            all_gc_samples.extend(chunk['gc_sample'])
        if 'quality_sample' in chunk:
            all_quality_samples.extend(chunk['quality_sample'])
    
    # Length statistics (should be uniform)
    all_lengths = [chunk['stats']['sequence_length']['mean'] for chunk in chunks]
    
    # Quality ranges
    min_quality = min(chunk['stats']['quality_scores']['min'] for chunk in chunks)
    max_quality = max(chunk['stats']['quality_scores']['max'] for chunk in chunks)
    
    file_stats = {
        'total_sequences': total_sequences,
        'num_chunks': len(chunks),
        'sequence_length': {
            'mean': np.mean(all_lengths) if all_lengths else 0,
            'min': min(all_lengths) if all_lengths else 0,
            'max': max(all_lengths) if all_lengths else 0
        },
        'gc_content': {
            'mean': mean_gc,
            'std': np.std(all_gc_samples) if all_gc_samples else 0,
            'samples': all_gc_samples[:10000]  # Limit for memory
        },
        'quality_scores': {
            'mean': mean_quality,
            'min': min_quality,
            'max': max_quality,
            'samples': all_quality_samples[:10000]  # Limit for memory
        },
        'chunk_summary': [
            {
                'chunk_number': i+1,
                'sequences': chunk['stats']['total_sequences'],
                'gc_content': chunk['stats']['gc_content']['mean'],
                'quality': chunk['stats']['quality_scores']['mean']
            }
            for i, chunk in enumerate(chunks)
        ]
    }
    
    return file_stats

def create_comprehensive_plots(file_name, file_stats, output_dir):
    """Create comprehensive FastQC-style plots for a file"""
    
    print(f"Creating plots for {file_name}...")
    
    # Set up the plot style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create a large figure with subplots
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle(f'Quality Analysis Report: {file_name}', fontsize=16, fontweight='bold')
    
    # 1. File Overview (top left)
    ax1 = plt.subplot(2, 3, 1)
    overview_text = f"""
File: {file_name}
Total Sequences: {file_stats['total_sequences']:,}
Chunks Processed: {file_stats['num_chunks']}
Sequence Length: {file_stats['sequence_length']['mean']:.0f} bp
Mean GC Content: {file_stats['gc_content']['mean']:.1f}%
Mean Quality: Q{file_stats['quality_scores']['mean']:.1f}
Quality Range: Q{file_stats['quality_scores']['min']}-Q{file_stats['quality_scores']['max']}
"""
    ax1.text(0.1, 0.9, overview_text, transform=ax1.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.axis('off')
    ax1.set_title('File Overview', fontweight='bold')
    
    # 2. GC Content Distribution (top middle)
    ax2 = plt.subplot(2, 3, 2)
    if file_stats['gc_content']['samples']:
        gc_data = file_stats['gc_content']['samples']
        ax2.hist(gc_data, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        ax2.axvline(file_stats['gc_content']['mean'], color='red', linestyle='--', linewidth=2,
                   label=f"Mean: {file_stats['gc_content']['mean']:.1f}%")
        
        # Add theoretical normal range
        ax2.axvspan(40, 60, alpha=0.1, color='green', label='Typical Range (40-60%)')
        
        ax2.set_xlabel('GC Content (%)')
        ax2.set_ylabel('Number of Reads')
        ax2.set_title('GC Content Distribution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    # 3. Quality Score Distribution (top right)
    ax3 = plt.subplot(2, 3, 3)
    if file_stats['quality_scores']['samples']:
        qual_data = file_stats['quality_scores']['samples']
        ax3.hist(qual_data, bins=range(0, 42), alpha=0.7, color='lightcoral', edgecolor='black')
        ax3.axvline(file_stats['quality_scores']['mean'], color='darkred', linestyle='--', linewidth=2,
                   label=f"Mean: Q{file_stats['quality_scores']['mean']:.1f}")
        
        # Add quality thresholds
        ax3.axvline(30, color='green', linestyle=':', label='Q30 (99.9% accuracy)')
        ax3.axvline(20, color='orange', linestyle=':', label='Q20 (99% accuracy)')
        ax3.axvline(10, color='red', linestyle=':', label='Q10 (90% accuracy)')
        
        ax3.set_xlabel('Quality Score (Phred)')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Quality Score Distribution')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # 4. Chunk-wise Quality Trends (bottom left)
    ax4 = plt.subplot(2, 3, 4)
    chunk_numbers = [chunk['chunk_number'] for chunk in file_stats['chunk_summary']]
    chunk_qualities = [chunk['quality'] for chunk in file_stats['chunk_summary']]
    
    ax4.plot(chunk_numbers, chunk_qualities, 'o-', color='purple', linewidth=2, markersize=4)
    ax4.axhline(30, color='green', linestyle='--', alpha=0.7, label='Q30')
    ax4.axhline(20, color='orange', linestyle='--', alpha=0.7, label='Q20')
    
    ax4.set_xlabel('Chunk Number')
    ax4.set_ylabel('Mean Quality Score')
    ax4.set_title('Quality Across File Chunks')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. Chunk-wise GC Content Trends (bottom middle)
    ax5 = plt.subplot(2, 3, 5)
    chunk_gc = [chunk['gc_content'] for chunk in file_stats['chunk_summary']]
    
    ax5.plot(chunk_numbers, chunk_gc, 's-', color='teal', linewidth=2, markersize=4)
    ax5.axhline(50, color='gray', linestyle='--', alpha=0.7, label='50% GC')
    ax5.axhspan(40, 60, alpha=0.1, color='green', label='Normal Range')
    
    ax5.set_xlabel('Chunk Number')
    ax5.set_ylabel('Mean GC Content (%)')
    ax5.set_title('GC Content Across File Chunks')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # 6. Quality Assessment Summary (bottom right)
    ax6 = plt.subplot(2, 3, 6)
    
    # Quality assessment
    assessment = []
    if file_stats['quality_scores']['mean'] >= 30:
        assessment.append("✅ Excellent Quality (Q≥30)")
    elif file_stats['quality_scores']['mean'] >= 20:
        assessment.append("✅ Good Quality (Q≥20)")
    else:
        assessment.append("⚠️ Low Quality (Q<20)")
    
    if 35 <= file_stats['gc_content']['mean'] <= 65:
        assessment.append("✅ Normal GC Content")
    else:
        assessment.append("⚠️ Unusual GC Content")
    
    if file_stats['sequence_length']['min'] == file_stats['sequence_length']['max']:
        assessment.append("✅ Uniform Read Length")
    else:
        assessment.append("⚠️ Variable Read Length")
    
    # Calculate some additional stats
    high_qual_percentage = sum(1 for q in file_stats['quality_scores']['samples'][:1000] if q >= 30) / 10
    assessment.append(f"📊 Q30 Bases: ~{high_qual_percentage:.1f}%")
    
    assessment_text = '\n'.join(assessment)
    
    ax6.text(0.1, 0.9, assessment_text, transform=ax6.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))
    ax6.set_xlim(0, 1)
    ax6.set_ylim(0, 1)
    ax6.axis('off')
    ax6.set_title('Quality Assessment', fontweight='bold')
    
    plt.tight_layout()
    
    # Save the plot
    safe_filename = file_name.replace('.fastq', '').replace('.', '_')
    output_file = output_dir / f"{safe_filename}_comprehensive_analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved: {output_file}")
    
    plt.show()

def generate_file_wise_report(file_results, output_dir):
    """Generate comprehensive file-wise report with statistics and plots"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print(f"\n=== Generating File-wise Analysis ===")
    
    # Process each file
    file_summaries = {}
    
    for file_name, chunks in file_results.items():
        print(f"\nProcessing {file_name}...")
        
        # Combine statistics from all chunks
        file_stats = combine_file_statistics(chunks)
        
        if file_stats:
            # Print summary
            print(f"  Total sequences: {file_stats['total_sequences']:,}")
            print(f"  Chunks: {file_stats['num_chunks']}")
            print(f"  Mean length: {file_stats['sequence_length']['mean']:.1f} bp")
            print(f"  Mean GC: {file_stats['gc_content']['mean']:.1f}%")
            print(f"  Mean quality: Q{file_stats['quality_scores']['mean']:.1f}")
            
            # Create comprehensive plots
            create_comprehensive_plots(file_name, file_stats, output_dir)
            
            # Save detailed file statistics
            stats_file = output_dir / f"{file_name.replace('.fastq', '')}_detailed_stats.json"
            with open(stats_file, 'w') as f:
                # Remove large sample arrays for JSON storage
                save_stats = file_stats.copy()
                save_stats['gc_content']['samples'] = save_stats['gc_content']['samples'][:1000]
                save_stats['quality_scores']['samples'] = save_stats['quality_scores']['samples'][:1000]
                json.dump(save_stats, f, indent=2)
            
            file_summaries[file_name] = file_stats
    
    # Create combined summary
    total_sequences = sum(stats['total_sequences'] for stats in file_summaries.values())
    total_chunks = sum(stats['num_chunks'] for stats in file_summaries.values())
    
    combined_summary = {
        'analysis_date': str(Path().cwd()),
        'total_files_analyzed': len(file_summaries),
        'total_chunks_processed': total_chunks,
        'total_sequences_analyzed': total_sequences,
        'file_summaries': {
            file_name: {
                'total_sequences': stats['total_sequences'],
                'mean_gc_content': stats['gc_content']['mean'],
                'mean_quality': stats['quality_scores']['mean'],
                'num_chunks': stats['num_chunks']
            }
            for file_name, stats in file_summaries.items()
        }
    }
    
    # Save combined summary
    summary_file = output_dir / "file_wise_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(combined_summary, f, indent=2)
    
    # Create text report
    report_file = output_dir / "file_wise_analysis_report.txt"
    with open(report_file, 'w') as f:
        f.write("FASTQ File-wise Analysis Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Analysis Summary:\n")
        f.write(f"  Files analyzed: {len(file_summaries)}\n")
        f.write(f"  Total chunks processed: {total_chunks}\n")
        f.write(f"  Total sequences: {total_sequences:,}\n\n")
        
        for file_name, stats in file_summaries.items():
            f.write(f"{file_name}:\n")
            f.write(f"  Sequences: {stats['total_sequences']:,}\n")
            f.write(f"  Chunks: {stats['num_chunks']}\n")
            f.write(f"  Mean length: {stats['sequence_length']['mean']:.1f} bp\n")
            f.write(f"  GC content: {stats['gc_content']['mean']:.1f}% ± {stats['gc_content']['std']:.1f}%\n")
            f.write(f"  Quality: Q{stats['quality_scores']['mean']:.1f} (range: Q{stats['quality_scores']['min']}-Q{stats['quality_scores']['max']})\n")
            f.write("\n")
    
    print(f"\n=== File-wise Analysis Complete ===")
    print(f"Results saved to: {output_dir}")
    print(f"Summary report: {report_file}")
    print(f"Detailed plots created for each file")
    
    return file_summaries

def main():
    """Main function"""
    
    # Set paths
    base_dir = Path(r"c:\Users\PAVILION\Desktop\genome")
    results_dir = base_dir / "chunk_results"
    output_dir = base_dir / "file_analysis"
    
    print("=== FASTQ File-wise Results Analysis ===")
    
    # Load chunk results
    file_results = load_chunk_results(results_dir)
    
    if not file_results:
        print("No chunk results found!")
        print("Make sure you've run the batch processing first.")
        return
    
    # Generate comprehensive file-wise analysis
    file_summaries = generate_file_wise_report(file_results, output_dir)
    
    print(f"\n=== Analysis Complete ===")
    print("Check the generated plots and reports in the file_analysis directory!")

if __name__ == "__main__":
    main()