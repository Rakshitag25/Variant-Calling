#!/usr/bin/env python3
"""
Compare original and trimmed reads to show the differences
"""

import os
from pathlib import Path
import random

def parse_fastq_entries(file_path, num_entries=10):
    """Parse FASTQ entries from file"""
    entries = []
    
    with open(file_path, 'r', encoding='utf-8') as f:
        while len(entries) < num_entries:
            try:
                header = f.readline().strip()
                if not header:
                    break
                    
                sequence = f.readline().strip()
                plus = f.readline().strip()
                quality = f.readline().strip()
                
                if header.startswith('@') and plus.startswith('+'):
                    entries.append({
                        'header': header,
                        'sequence': sequence,
                        'quality': quality,
                        'length': len(sequence)
                    })
            except:
                break
    
    return entries

def get_quality_stats(quality_string):
    """Get quality statistics from quality string"""
    if not quality_string:
        return {'min': 0, 'max': 0, 'avg': 0}
    
    # Convert ASCII to Phred scores (assuming Phred+33)
    phred_scores = [ord(char) - 33 for char in quality_string]
    
    return {
        'min': min(phred_scores),
        'max': max(phred_scores),
        'avg': sum(phred_scores) / len(phred_scores)
    }

def find_differences(original_seq, trimmed_seq):
    """Find what was trimmed from the sequence"""
    if original_seq == trimmed_seq:
        return "No changes"
    
    # Check if it's a simple trimming from ends
    if trimmed_seq in original_seq:
        start_pos = original_seq.find(trimmed_seq)
        end_pos = start_pos + len(trimmed_seq)
        
        trimmed_left = original_seq[:start_pos] if start_pos > 0 else ""
        trimmed_right = original_seq[end_pos:] if end_pos < len(original_seq) else ""
        
        changes = []
        if trimmed_left:
            changes.append(f"Left: '{trimmed_left}' ({len(trimmed_left)} bp)")
        if trimmed_right:
            changes.append(f"Right: '{trimmed_right}' ({len(trimmed_right)} bp)")
        
        return "Trimmed - " + ", ".join(changes) if changes else "No changes"
    else:
        return f"Modified (original: {len(original_seq)} bp → trimmed: {len(trimmed_seq)} bp)"

def compare_reads():
    """Compare original and trimmed reads"""
    
    # File paths
    original_chunk = "ERR2304551_1.fastq/chunks/ERR2304551_1_chunk_001.fastq"
    trimmed_file = "fastpandtrim/ERR2304551_1_trimmed.fastq"
    
    if not os.path.exists(original_chunk):
        print(f"❌ Original chunk file not found: {original_chunk}")
        return
    
    if not os.path.exists(trimmed_file):
        print(f"❌ Trimmed file not found: {trimmed_file}")
        return
    
    print("=" * 80)
    print("ORIGINAL vs TRIMMED READS COMPARISON")
    print("=" * 80)
    
    # Get sample reads from original chunk
    print(f"📖 Reading sample reads from: {original_chunk}")
    original_reads = parse_fastq_entries(original_chunk, 20)
    
    print(f"📖 Reading sample reads from: {trimmed_file}")
    trimmed_reads = parse_fastq_entries(trimmed_file, 20)
    
    if not original_reads:
        print("❌ No reads found in original file")
        return
    
    if not trimmed_reads:
        print("❌ No reads found in trimmed file")
        return
    
    print(f"✅ Found {len(original_reads)} original reads and {len(trimmed_reads)} trimmed reads")
    print()
    
    # Compare first few reads that should correspond
    for i, (orig, trimmed) in enumerate(zip(original_reads[:5], trimmed_reads[:5])):
        print(f"🔍 READ {i+1} COMPARISON")
        print("-" * 60)
        
        # Extract read ID (remove @)
        orig_id = orig['header'][1:].split()[0]
        trimmed_id = trimmed['header'][1:].split()[0]
        
        print(f"Original ID:  {orig_id}")
        print(f"Trimmed ID:   {trimmed_id}")
        print(f"IDs match:    {'✅ Yes' if orig_id == trimmed_id else '❌ No'}")
        print()
        
        # Length comparison
        print(f"Length:       {orig['length']} bp → {trimmed['length']} bp (change: {trimmed['length'] - orig['length']:+d})")
        
        # Quality statistics
        orig_qual = get_quality_stats(orig['quality'])
        trimmed_qual = get_quality_stats(trimmed['quality'])
        
        print(f"Quality (avg): Q{orig_qual['avg']:.1f} → Q{trimmed_qual['avg']:.1f}")
        print(f"Quality (min): Q{orig_qual['min']} → Q{trimmed_qual['min']}")
        print(f"Quality (max): Q{orig_qual['max']} → Q{trimmed_qual['max']}")
        print()
        
        # Show sequences
        print("ORIGINAL SEQUENCE:")
        print(f"  {orig['sequence'][:80]}{'...' if len(orig['sequence']) > 80 else ''}")
        print()
        print("TRIMMED SEQUENCE:")
        print(f"  {trimmed['sequence'][:80]}{'...' if len(trimmed['sequence']) > 80 else ''}")
        print()
        
        # Show quality scores
        print("ORIGINAL QUALITY:")
        print(f"  {orig['quality'][:80]}{'...' if len(orig['quality']) > 80 else ''}")
        print()
        print("TRIMMED QUALITY:")
        print(f"  {trimmed['quality'][:80]}{'...' if len(trimmed['quality']) > 80 else ''}")
        print()
        
        # Analyze changes
        changes = find_differences(orig['sequence'], trimmed['sequence'])
        print(f"CHANGES: {changes}")
        
        # Show end regions if trimmed
        if orig['length'] != trimmed['length']:
            print("\nEND REGIONS ANALYSIS:")
            print(f"Original 5': {orig['sequence'][:10]}... (Q: {orig['quality'][:10]}...)")
            print(f"Original 3': ...{orig['sequence'][-10:]} (Q: ...{orig['quality'][-10:]})")
            print(f"Trimmed 5':  {trimmed['sequence'][:10]}... (Q: {trimmed['quality'][:10]}...)")
            print(f"Trimmed 3':  ...{trimmed['sequence'][-10:]} (Q: ...{trimmed['quality'][-10:]})")
        
        print("=" * 60)
        print()
    
    # Summary statistics
    print("📊 SUMMARY STATISTICS")
    print("-" * 40)
    
    orig_lengths = [read['length'] for read in original_reads]
    trimmed_lengths = [read['length'] for read in trimmed_reads]
    
    orig_qualities = [get_quality_stats(read['quality'])['avg'] for read in original_reads]
    trimmed_qualities = [get_quality_stats(read['quality'])['avg'] for read in trimmed_reads]
    
    print(f"Sample size:           {len(original_reads)} reads")
    print(f"Average length:        {sum(orig_lengths)/len(orig_lengths):.1f} → {sum(trimmed_lengths)/len(trimmed_lengths):.1f} bp")
    print(f"Length range:          {min(orig_lengths)}-{max(orig_lengths)} → {min(trimmed_lengths)}-{max(trimmed_lengths)} bp")
    print(f"Average quality:       Q{sum(orig_qualities)/len(orig_qualities):.1f} → Q{sum(trimmed_qualities)/len(trimmed_qualities):.1f}")
    
    # Count changes
    unchanged = sum(1 for orig, trimmed in zip(original_reads, trimmed_reads) if orig['sequence'] == trimmed['sequence'])
    changed = len(original_reads) - unchanged
    
    print(f"Unchanged reads:       {unchanged}/{len(original_reads)} ({unchanged/len(original_reads)*100:.1f}%)")
    print(f"Modified reads:        {changed}/{len(original_reads)} ({changed/len(original_reads)*100:.1f}%)")

def main():
    print("Starting read comparison analysis...")
    compare_reads()

if __name__ == "__main__":
    main()