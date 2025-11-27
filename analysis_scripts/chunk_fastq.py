#!/usr/bin/env python3
"""
Chunk-based FASTQ Processing
Breaks large files into manageable chunks for processing
"""

import os
from pathlib import Path
import subprocess

def split_fastq_into_chunks(input_file, chunk_size=1000000, output_dir=None):
    """
    Split large FASTQ file into smaller chunks
    chunk_size: number of reads per chunk
    """
    
    if output_dir is None:
        output_dir = Path(input_file).parent / "chunks"
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print(f"Splitting {input_file} into chunks of {chunk_size:,} reads...")
    
    chunk_files = []
    current_chunk = 1
    current_read_count = 0
    current_chunk_file = None
    
    try:
        with open(input_file, 'r') as f:
            line_count = 0
            
            for line in f:
                line_count += 1
                position_in_read = (line_count - 1) % 4
                
                # Start new chunk if needed
                if position_in_read == 0 and current_read_count % chunk_size == 0:
                    if current_chunk_file:
                        current_chunk_file.close()
                    
                    chunk_filename = output_dir / f"{Path(input_file).stem}_chunk_{current_chunk:03d}.fastq"
                    current_chunk_file = open(chunk_filename, 'w')
                    chunk_files.append(chunk_filename)
                    print(f"Creating chunk {current_chunk}: {chunk_filename}")
                    current_chunk += 1
                
                current_chunk_file.write(line)
                
                # Count completed reads
                if position_in_read == 3:
                    current_read_count += 1
                    
                    if current_read_count % 100000 == 0:
                        print(f"Processed {current_read_count:,} reads...")
        
        if current_chunk_file:
            current_chunk_file.close()
            
    except Exception as e:
        print(f"Error splitting file: {e}")
        return []
    
    print(f"Split complete: {len(chunk_files)} chunks created")
    return chunk_files

def process_chunks_sequentially(chunk_files, analysis_script):
    """Process each chunk sequentially and combine results"""
    
    all_results = []
    
    for i, chunk_file in enumerate(chunk_files):
        print(f"\nProcessing chunk {i+1}/{len(chunk_files)}: {chunk_file}")
        
        try:
            # Run analysis on chunk
            result = subprocess.run([
                "python", analysis_script, str(chunk_file)
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"✅ Chunk {i+1} processed successfully")
                all_results.append(result.stdout)
            else:
                print(f"❌ Error processing chunk {i+1}: {result.stderr}")
                
        except Exception as e:
            print(f"❌ Exception processing chunk {i+1}: {e}")
    
    return all_results

def estimate_optimal_chunk_size(file_path, target_memory_mb=500):
    """Estimate optimal chunk size based on available memory"""
    
    file_size = os.path.getsize(file_path)
    file_size_mb = file_size / (1024 * 1024)
    
    print(f"File size: {file_size_mb:.1f} MB")
    
    # Estimate reads per MB (roughly 2500 reads per MB for typical FASTQ)
    estimated_reads_per_mb = 2500
    total_estimated_reads = int(file_size_mb * estimated_reads_per_mb)
    
    # Calculate chunk size to stay within target memory
    target_memory_ratio = target_memory_mb / file_size_mb
    chunk_size = max(10000, int(total_estimated_reads * target_memory_ratio))
    
    num_chunks = (total_estimated_reads + chunk_size - 1) // chunk_size
    
    print(f"Estimated total reads: {total_estimated_reads:,}")
    print(f"Recommended chunk size: {chunk_size:,} reads")
    print(f"Estimated number of chunks: {num_chunks}")
    
    return chunk_size

def main():
    """Main chunking workflow"""
    
    # Configuration
    large_files = [
        Path(r"c:\Users\PAVILION\Desktop\genome\ERR2304551_1.fastq\ERR2304551_1.fastq"),
        Path(r"c:\Users\PAVILION\Desktop\genome\ERR2304551_2.fastq\ERR2304551_2.fastq"),
    ]
    
    for file_path in large_files:
        if not file_path.exists():
            print(f"File not found: {file_path}")
            continue
            
        print(f"\n{'='*60}")
        print(f"Processing: {file_path}")
        
        # Estimate optimal chunk size
        chunk_size = estimate_optimal_chunk_size(file_path)
        
        # Ask user for confirmation
        user_input = input(f"\nProceed with chunking? (y/n): ")
        if user_input.lower() != 'y':
            print("Skipping file...")
            continue
        
        # Split file into chunks
        chunk_files = split_fastq_into_chunks(file_path, chunk_size)
        
        if chunk_files:
            print(f"\nChunking complete! Created {len(chunk_files)} chunks")
            print("You can now process each chunk individually with your analysis script.")
            print("\nChunk files created:")
            for chunk_file in chunk_files:
                file_size = os.path.getsize(chunk_file) / (1024 * 1024)
                print(f"  {chunk_file.name} ({file_size:.1f} MB)")

if __name__ == "__main__":
    main()