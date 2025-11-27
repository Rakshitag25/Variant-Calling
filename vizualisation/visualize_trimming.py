#!/usr/bin/env python3
"""
Generate visualization plots for the trimming results
"""

import json
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path
import pandas as pd

def load_trimming_data():
    """Load the comprehensive trimming report"""
    report_file = Path("fastpandtrim") / "ERR2304551_1_final_trim_report.json"
    
    if not report_file.exists():
        print(f"Report file not found: {report_file}")
        return None
    
    with open(report_file, 'r') as f:
        data = json.load(f)
    
    return data

def create_trimming_plots(data):
    """Create comprehensive plots for trimming results"""
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 16))
    
    # 1. Overall Summary (Top Left)
    ax1 = plt.subplot(3, 3, 1)
    summary = data['summary']
    labels = ['Passed', 'Too Short', 'Low Quality']
    sizes = [summary['passed_reads'], summary['too_short'], summary['low_quality']]
    colors = ['#2ecc71', '#f39c12', '#e74c3c']
    
    wedges, texts, autotexts = ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.2f%%', startangle=90)
    ax1.set_title(f'Read Filtering Results\nTotal: {summary["total_reads"]:,} reads', fontsize=12, fontweight='bold')
    
    # 2. Length Distribution (Top Middle)
    ax2 = plt.subplot(3, 3, 2)
    length_stats = data['length_statistics']
    categories = ['Before Trimming', 'After Trimming']
    lengths = [length_stats['avg_length_before'], length_stats['avg_length_after']]
    colors = ['#3498db', '#27ae60']
    
    bars = ax2.bar(categories, lengths, color=colors, alpha=0.8)
    ax2.set_ylabel('Average Read Length (bp)')
    ax2.set_title('Average Read Length\nBefore vs After Trimming', fontsize=12, fontweight='bold')
    ax2.set_ylim(148, 151)
    
    # Add value labels on bars
    for bar, length in zip(bars, lengths):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05, 
                f'{length:.1f} bp', ha='center', va='bottom', fontweight='bold')
    
    # 3. Processing Statistics (Top Right)
    ax3 = plt.subplot(3, 3, 3)
    proc_stats = data['processing']
    metrics = ['Chunks\nProcessed', 'Processing\nTime (hrs)', 'Processing\nRate (K reads/s)']
    values = [proc_stats['chunks_processed'], 
              proc_stats['total_processing_time_seconds']/3600,
              proc_stats['total_processing_time_seconds']/data['summary']['total_reads']*1000]
    
    bars = ax3.bar(metrics, values, color=['#9b59b6', '#e67e22', '#1abc9c'])
    ax3.set_title('Processing Metrics', fontsize=12, fontweight='bold')
    
    # Add value labels
    for bar, value in zip(bars, values):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(values)*0.02, 
                f'{value:.1f}', ha='center', va='bottom', fontweight='bold')
    
    # 4. Bases Processed (Middle Left)
    ax4 = plt.subplot(3, 3, 4)
    bases_before = length_stats['input_bases'] / 1e9  # Convert to Gb
    bases_after = length_stats['output_bases'] / 1e9
    bases_removed = length_stats['bases_removed'] / 1e6  # Convert to Mb
    
    categories = ['Input\n(Gb)', 'Output\n(Gb)', 'Removed\n(Mb)']
    values = [bases_before, bases_after, bases_removed]
    colors = ['#3498db', '#27ae60', '#e74c3c']
    
    bars = ax4.bar(categories, values, color=colors, alpha=0.8)
    ax4.set_title('Base Processing Statistics', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Bases')
    
    # Add value labels
    for bar, value, unit in zip(bars, values, ['Gb', 'Gb', 'Mb']):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(values)*0.02, 
                f'{value:.1f} {unit}', ha='center', va='bottom', fontweight='bold')
    
    # 5. Per-Chunk Performance (Middle)
    ax5 = plt.subplot(3, 3, (5, 6))  # Span two columns
    chunk_details = pd.DataFrame(data['chunk_details'])
    chunk_details['chunk_num'] = range(1, len(chunk_details) + 1)
    
    # Plot pass rate per chunk
    ax5.plot(chunk_details['chunk_num'], chunk_details['pass_rate'], 
            marker='o', linewidth=2, markersize=4, color='#2ecc71')
    ax5.set_xlabel('Chunk Number')
    ax5.set_ylabel('Pass Rate (%)')
    ax5.set_title('Read Pass Rate Per Chunk', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    ax5.set_ylim(99.5, 100.1)
    
    # 6. Quality Distribution (Bottom Left)
    ax6 = plt.subplot(3, 3, 7)
    filter_reasons = ['Too Short', 'Low Quality']
    filter_counts = [summary['too_short'], summary['low_quality']]
    colors = ['#f39c12', '#e74c3c']
    
    bars = ax6.bar(filter_reasons, filter_counts, color=colors, alpha=0.8)
    ax6.set_ylabel('Number of Reads')
    ax6.set_title('Filtering Breakdown', fontsize=12, fontweight='bold')
    
    # Add value labels
    for bar, count in zip(bars, filter_counts):
        ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(filter_counts)*0.02, 
                f'{count:,}', ha='center', va='bottom', fontweight='bold')
    
    # 7. Length Change Distribution (Bottom Middle)
    ax7 = plt.subplot(3, 3, 8)
    chunk_details['length_change'] = chunk_details['avg_length_before'] - chunk_details['avg_length_after']
    
    ax7.hist(chunk_details['length_change'], bins=15, color='#3498db', alpha=0.7, edgecolor='black')
    ax7.set_xlabel('Average Length Change (bp)')
    ax7.set_ylabel('Number of Chunks')
    ax7.set_title('Distribution of Length Changes', fontsize=12, fontweight='bold')
    ax7.grid(True, alpha=0.3)
    
    # 8. Summary Statistics (Bottom Right)
    ax8 = plt.subplot(3, 3, 9)
    ax8.axis('off')
    
    # Create summary text
    summary_text = f"""
TRIMMING SUMMARY - ERR2304551_1

Total Reads: {summary['total_reads']:,}
Passed: {summary['passed_reads']:,} ({summary['pass_rate']:.1f}%)
Filtered: {summary['filtered_out']:,} ({100-summary['pass_rate']:.1f}%)

Average Length:
  Before: {length_stats['avg_length_before']:.1f} bp
  After: {length_stats['avg_length_after']:.1f} bp
  Trimmed: {length_stats['length_reduction']:.1f} bp

Processing:
  Time: {proc_stats['total_processing_time_minutes']/60:.1f} hours
  Rate: {data['summary']['total_reads']/proc_stats['total_processing_time_seconds']:.0f} reads/sec
  
Quality Control:
  Too Short: {summary['too_short']:,} reads
  Low Quality: {summary['low_quality']:,} reads
  
File Size: ~14.0 GB trimmed
"""
    
    ax8.text(0.05, 0.95, summary_text, transform=ax8.transAxes, fontsize=11,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8))
    
    plt.tight_layout(pad=3.0)
    
    # Save the plot
    output_file = Path("fastpandtrim") / "ERR2304551_1_trimming_summary.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"✓ Trimming summary plot saved: {output_file}")

def main():
    print("=== Generating Trimming Visualization ===")
    
    # Load data
    data = load_trimming_data()
    if data is None:
        return
    
    print("✓ Loaded trimming report data")
    
    # Create plots
    create_trimming_plots(data)
    
    print("✓ Trimming visualization completed")

if __name__ == "__main__":
    main()