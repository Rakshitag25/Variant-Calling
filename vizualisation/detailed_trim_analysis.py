#!/usr/bin/env python3
"""
Detailed analysis of specific trimming examples showing before/after with quality visualization
"""

import os
from pathlib import Path

def get_quality_char_meaning(phred_score):
    """Convert Phred score to quality meaning"""
    if phred_score >= 30:
        return "Excellent"
    elif phred_score >= 20:
        return "Good"
    elif phred_score >= 10:
        return "Poor"
    else:
        return "Bad"

def visualize_quality_scores(quality_string, sequence, title):
    """Create a visual representation of quality scores"""
    print(f"\n{title}")
    print("=" * len(title))
    
    # Convert ASCII to Phred scores
    phred_scores = [ord(char) - 33 for char in quality_string]
    
    # Create quality visualization
    print("Position: ", end="")
    for i in range(0, len(sequence), 10):
        print(f"{i:>10}", end="")
    print()
    
    print("Sequence: ", end="")
    for i, base in enumerate(sequence):
        if i % 10 == 0 and i > 0:
            print(f"{base}", end="")
        else:
            print(f"{base}", end="")
    print()
    
    print("Quality:  ", end="")
    for i, (char, score) in enumerate(zip(quality_string, phred_scores)):
        if score >= 30:
            symbol = "+"  # Excellent
        elif score >= 20:
            symbol = "="  # Good
        elif score >= 10:
            symbol = "-"  # Poor
        else:
            symbol = "x"  # Bad
        print(symbol, end="")
    print()
    
    print("Phred Q:  ", end="")
    for i, score in enumerate(phred_scores):
        if i % 5 == 0:
            print(f"{score:>2}", end="")
        else:
            print("  ", end="")
    print()
    
    # Quality legend
    print("\nQuality Legend: + Excellent (Q≥30) | = Good (Q≥20) | - Poor (Q≥10) | x Bad (Q<10)")
    
    # Quality statistics
    min_q = min(phred_scores)
    max_q = max(phred_scores)
    avg_q = sum(phred_scores) / len(phred_scores)
    
    print(f"Quality Stats: Min=Q{min_q} | Max=Q{max_q} | Avg=Q{avg_q:.1f}")
    
    return phred_scores

def analyze_specific_reads():
    """Analyze specific reads that show different types of trimming"""
    
    original_chunk = "ERR2304551_1.fastq/chunks/ERR2304551_1_chunk_001.fastq"
    trimmed_file = "fastpandtrim/ERR2304551_1_trimmed.fastq"
    
    print("🔬 DETAILED TRIMMING ANALYSIS")
    print("=" * 80)
    
    # Read first 100 reads to find interesting examples
    original_reads = []
    with open(original_chunk, 'r', encoding='utf-8') as f:
        for i in range(100):
            try:
                header = f.readline().strip()
                sequence = f.readline().strip()
                plus = f.readline().strip()
                quality = f.readline().strip()
                
                if header.startswith('@') and plus.startswith('+'):
                    original_reads.append({
                        'header': header,
                        'sequence': sequence,
                        'quality': quality
                    })
            except:
                break
    
    trimmed_reads = []
    with open(trimmed_file, 'r', encoding='utf-8') as f:
        for i in range(100):
            try:
                header = f.readline().strip()
                sequence = f.readline().strip()
                plus = f.readline().strip()
                quality = f.readline().strip()
                
                if header.startswith('@') and plus.startswith('+'):
                    trimmed_reads.append({
                        'header': header,
                        'sequence': sequence,
                        'quality': quality
                    })
            except:
                break
    
    # Find examples of different types of trimming
    examples = {
        'no_change': [],
        'minor_trim': [],
        'major_trim': [],
        'quality_trim': []
    }
    
    for i, (orig, trimmed) in enumerate(zip(original_reads, trimmed_reads)):
        length_diff = len(orig['sequence']) - len(trimmed['sequence'])
        
        if length_diff == 0:
            examples['no_change'].append((i, orig, trimmed))
        elif 1 <= length_diff <= 3:
            examples['minor_trim'].append((i, orig, trimmed))
        elif length_diff > 10:
            examples['major_trim'].append((i, orig, trimmed))
        elif length_diff > 3:
            examples['quality_trim'].append((i, orig, trimmed))
    
    # Show examples
    print("📊 TRIMMING CATEGORIES FOUND:")
    print(f"  No changes:    {len(examples['no_change'])} reads")
    print(f"  Minor trimming (1-3 bp): {len(examples['minor_trim'])} reads")
    print(f"  Quality trimming (4-10 bp): {len(examples['quality_trim'])} reads")
    print(f"  Major trimming (>10 bp): {len(examples['major_trim'])} reads")
    print()
    
    # Example 1: No change
    if examples['no_change']:
        idx, orig, trimmed = examples['no_change'][0]
        print(f"🟢 EXAMPLE 1: HIGH QUALITY READ (No trimming needed)")
        print(f"Read #{idx+1} - {orig['header']}")
        print(f"Length: {len(orig['sequence'])} bp (unchanged)")
        
        visualize_quality_scores(orig['quality'], orig['sequence'], "Quality Profile")
        
        print(f"\n✅ This read maintained excellent quality throughout, so no trimming was needed.")
        print("-" * 80)
    
    # Example 2: Minor trimming
    if examples['minor_trim']:
        idx, orig, trimmed = examples['minor_trim'][0]
        print(f"\n🟡 EXAMPLE 2: MINOR TRIMMING")
        print(f"Read #{idx+1} - {orig['header']}")
        print(f"Length: {len(orig['sequence'])} bp → {len(trimmed['sequence'])} bp (trimmed {len(orig['sequence'])-len(trimmed['sequence'])} bp)")
        
        print(f"\nORIGINAL READ:")
        visualize_quality_scores(orig['quality'], orig['sequence'], "Original Quality Profile")
        
        print(f"\nTRIMMED READ:")
        visualize_quality_scores(trimmed['quality'], trimmed['sequence'], "Trimmed Quality Profile")
        
        # Show what was trimmed
        if len(orig['sequence']) > len(trimmed['sequence']):
            if trimmed['sequence'] in orig['sequence']:
                start = orig['sequence'].find(trimmed['sequence'])
                if start > 0:
                    print(f"\n❌ TRIMMED FROM 5' END: '{orig['sequence'][:start]}'")
                    print(f"   Quality scores: {orig['quality'][:start]} (Phred: {[ord(c)-33 for c in orig['quality'][:start]]})")
                
                end = start + len(trimmed['sequence'])
                if end < len(orig['sequence']):
                    print(f"\n❌ TRIMMED FROM 3' END: '{orig['sequence'][end:]}'")
                    print(f"   Quality scores: {orig['quality'][end:]} (Phred: {[ord(c)-33 for c in orig['quality'][end:]]})")
        
        print("-" * 80)
    
    # Example 3: Major trimming
    if examples['major_trim']:
        idx, orig, trimmed = examples['major_trim'][0]
        print(f"\n🔴 EXAMPLE 3: MAJOR TRIMMING (Low quality regions)")
        print(f"Read #{idx+1} - {orig['header']}")
        print(f"Length: {len(orig['sequence'])} bp → {len(trimmed['sequence'])} bp (trimmed {len(orig['sequence'])-len(trimmed['sequence'])} bp)")
        
        print(f"\nORIGINAL READ:")
        orig_scores = visualize_quality_scores(orig['quality'], orig['sequence'], "Original Quality Profile")
        
        print(f"\nTRIMMED READ:")
        trimmed_scores = visualize_quality_scores(trimmed['quality'], trimmed['sequence'], "Trimmed Quality Profile")
        
        # Analyze the trimming decision
        if len(orig['sequence']) > len(trimmed['sequence']):
            trimmed_length = len(orig['sequence']) - len(trimmed['sequence'])
            
            # Check quality at trimmed regions
            if trimmed['sequence'] in orig['sequence']:
                start = orig['sequence'].find(trimmed['sequence'])
                end = start + len(trimmed['sequence'])
                
                if start > 0:
                    trimmed_5_qual = [ord(c)-33 for c in orig['quality'][:start]]
                    print(f"\n❌ 5' END TRIMMED: '{orig['sequence'][:start]}'")
                    print(f"   Average quality: Q{sum(trimmed_5_qual)/len(trimmed_5_qual):.1f}")
                
                if end < len(orig['sequence']):
                    trimmed_3_qual = [ord(c)-33 for c in orig['quality'][end:]]
                    print(f"\n❌ 3' END TRIMMED: '{orig['sequence'][end:]}'")
                    print(f"   Average quality: Q{sum(trimmed_3_qual)/len(trimmed_3_qual):.1f}")
                
                kept_qual = [ord(c)-33 for c in trimmed['quality']]
                print(f"\n✅ KEPT REGION average quality: Q{sum(kept_qual)/len(kept_qual):.1f}")
        
        print("-" * 80)
    
    # Summary
    print(f"\n📈 OVERALL TRIMMING EFFECTIVENESS:")
    total_reads = len(original_reads)
    total_orig_length = sum(len(read['sequence']) for read in original_reads)
    total_trim_length = sum(len(read['sequence']) for read in trimmed_reads)
    
    print(f"  Sample reads analyzed: {total_reads}")
    print(f"  Total bases before: {total_orig_length:,} bp")
    print(f"  Total bases after: {total_trim_length:,} bp")
    print(f"  Bases removed: {total_orig_length - total_trim_length:,} bp ({(total_orig_length-total_trim_length)/total_orig_length*100:.2f}%)")
    print(f"  Data retention: {total_trim_length/total_orig_length*100:.2f}%")
    
    # Quality improvement
    orig_avg_qual = []
    trim_avg_qual = []
    
    for orig, trimmed in zip(original_reads, trimmed_reads):
        orig_scores = [ord(c)-33 for c in orig['quality']]
        trim_scores = [ord(c)-33 for c in trimmed['quality']]
        orig_avg_qual.append(sum(orig_scores)/len(orig_scores))
        trim_avg_qual.append(sum(trim_scores)/len(trim_scores))
    
    print(f"  Average quality before: Q{sum(orig_avg_qual)/len(orig_avg_qual):.1f}")
    print(f"  Average quality after: Q{sum(trim_avg_qual)/len(trim_avg_qual):.1f}")
    print(f"  Quality improvement: +{sum(trim_avg_qual)/len(trim_avg_qual) - sum(orig_avg_qual)/len(orig_avg_qual):.1f} Phred points")

def main():
    print("Starting detailed trimming analysis...")
    analyze_specific_reads()

if __name__ == "__main__":
    main()