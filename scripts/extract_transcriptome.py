#!/usr/bin/env python3
"""
Extract transcriptome sequences from genome FASTA and GTF annotation
用于从基因组 FASTA 和 GTF 注释文件中提取转录本序列

修复版：支持多种 GTF 格式，更健壮的错误处理
"""

import sys
import gzip
from collections import defaultdict

def parse_gtf(gtf_file):
    """
    解析 GTF 文件，提取外显子信息
    返回: {transcript_id: {'chr': chr, 'strand': strand, 'gene_name': name, 'exons': [(start, end), ...]}}
    """
    transcripts = defaultdict(lambda: {'chr': None, 'strand': None, 'gene_name': None, 'exons': []})
    
    print(f"Parsing GTF file: {gtf_file}")
    
    # 判断是否为压缩文件
    open_func = gzip.open if gtf_file.endswith('.gz') else open
    
    line_count = 0
    exon_count = 0
    feature_types = set()
    
    with open_func(gtf_file, 'rt') as f:
        for line in f:
            line_count += 1
            if line_count % 100000 == 0:
                print(f"  Processed {line_count} lines, {exon_count} exons...")
            
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            feature_types.add(feature_type)
            
            # 处理外显子和 CDS（有些 GTF 文件用 CDS 代替 exon）
            if feature_type not in ['exon', 'CDS']:
                continue
            
            chrom = fields[0]
            start = int(fields[3])  # GTF 是 1-based
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            # 解析属性字段（支持多种格式）
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                
                # 处理两种格式：key "value" 和 key=value
                if ' "' in attr:
                    key, value = attr.split(' "', 1)
                    attr_dict[key.strip()] = value.rstrip('"')
                elif '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key.strip()] = value.strip().strip('"')
            
            # 提取 transcript_id 和 gene_name
            transcript_id = attr_dict.get('transcript_id', '')
            gene_id = attr_dict.get('gene_id', '')
            gene_name = attr_dict.get('gene_name', gene_id)
            
            if not transcript_id:
                # 如果没有 transcript_id，尝试用 gene_id 作为 transcript_id
                transcript_id = gene_id
            
            if not transcript_id:
                continue
            
            # 存储外显子信息
            if transcripts[transcript_id]['chr'] is None:
                transcripts[transcript_id]['chr'] = chrom
                transcripts[transcript_id]['strand'] = strand
                transcripts[transcript_id]['gene_name'] = gene_name if gene_name else gene_id
            
            transcripts[transcript_id]['exons'].append((start, end))
            exon_count += 1
    
    print(f"  Total lines: {line_count}")
    print(f"  Feature types found: {', '.join(sorted(feature_types))}")
    print(f"  Total exons/CDS: {exon_count}")
    print(f"  Total transcripts: {len(transcripts)}")
    
    # 对每个转录本的外显子按位置排序并去重
    for transcript_id in transcripts:
        # 去重（有些 GTF 可能有重复的外显子）
        exons = list(set(transcripts[transcript_id]['exons']))
        exons.sort()
        transcripts[transcript_id]['exons'] = exons
    
    return dict(transcripts)


def load_genome(fasta_file):
    """
    加载基因组序列
    返回: {chr_id: sequence}
    """
    genome = {}
    current_chr = None
    current_seq = []
    
    print(f"\nLoading genome FASTA: {fasta_file}")
    
    # 判断是否为压缩文件
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    
    chr_count = 0
    with open_func(fasta_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 保存上一条染色体
                if current_chr is not None:
                    genome[current_chr] = ''.join(current_seq)
                    chr_count += 1
                    if chr_count <= 10:  # 只打印前 10 条
                        print(f"  Loaded {current_chr}: {len(genome[current_chr]):,} bp")
                
                # 开始新染色体
                # 处理多种 header 格式：">1 dna:chromosome..." 或 ">Chr1" 等
                header = line[1:].split()[0]
                current_chr = header
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        # 保存最后一条染色体
        if current_chr is not None:
            genome[current_chr] = ''.join(current_seq)
            chr_count += 1
            if chr_count <= 10:
                print(f"  Loaded {current_chr}: {len(genome[current_chr]):,} bp")
    
    if chr_count > 10:
        print(f"  ... and {chr_count - 10} more chromosomes")
    
    print(f"  Total chromosomes/scaffolds: {len(genome)}")
    return genome


def reverse_complement(seq):
    """反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def extract_transcript_sequence(genome, transcript_info):
    """
    从基因组中提取转录本序列
    """
    chrom = transcript_info['chr']
    strand = transcript_info['strand']
    exons = transcript_info['exons']
    
    if chrom not in genome:
        return None
    
    # 提取并拼接所有外显子序列
    transcript_seq = []
    for start, end in exons:
        # GTF 是 1-based, Python 是 0-based
        exon_seq = genome[chrom][start-1:end]
        transcript_seq.append(exon_seq)
    
    transcript_seq = ''.join(transcript_seq)
    
    # 如果是负链，取反向互补
    if strand == '-':
        transcript_seq = reverse_complement(transcript_seq)
    
    return transcript_seq


def write_transcriptome(transcripts, genome, output_file, gene_filter=None):
    """
    提取并写入转录组序列
    
    Args:
        transcripts: 转录本信息字典
        genome: 基因组序列字典
        output_file: 输出文件路径
        gene_filter: 可选，只提取特定基因（用于快速测试）
    """
    print(f"\nExtracting transcriptome sequences...")
    
    if gene_filter:
        print(f"  Gene filter: '{gene_filter}'")
    
    success_count = 0
    fail_count = 0
    fail_reasons = defaultdict(int)
    
    with open(output_file, 'w') as out:
        for transcript_id, info in transcripts.items():
            gene_name = info['gene_name']
            
            # 如果指定了基因过滤，只处理该基因
            if gene_filter:
                # 支持部分匹配（不区分大小写）
                if not (gene_filter.upper() in gene_name.upper() or 
                       gene_filter.upper() in transcript_id.upper()):
                    continue
            
            seq = extract_transcript_sequence(genome, info)
            
            if seq is None:
                fail_count += 1
                fail_reasons[f"Chromosome '{info['chr']}' not found"] += 1
                continue
            
            if len(seq) == 0:
                fail_count += 1
                fail_reasons["Empty sequence"] += 1
                continue
            
            # 写入 FASTA 格式
            # 格式: >transcript_id gene_name chr:start-end(strand)
            if info['exons']:
                location = f"{info['chr']}:{info['exons'][0][0]}-{info['exons'][-1][1]}({info['strand']})"
            else:
                location = f"{info['chr']}:unknown"
            
            header = f">{transcript_id} {gene_name} {location}"
            out.write(header + '\n')
            
            # 每行 60 个字符
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + '\n')
            
            success_count += 1
            
            if success_count % 1000 == 0:
                print(f"  Extracted {success_count} transcripts...")
    
    print(f"\n✓ Successfully extracted {success_count} transcripts")
    
    if fail_count > 0:
        print(f"✗ Failed to extract {fail_count} transcripts:")
        for reason, count in fail_reasons.items():
            print(f"    - {reason}: {count}")
    
    print(f"\nOutput saved to: {output_file}")
    
    # 显示文件大小
    import os
    file_size = os.path.getsize(output_file)
    if file_size > 1024*1024:
        print(f"File size: {file_size / (1024*1024):.1f} MB")
    else:
        print(f"File size: {file_size / 1024:.1f} KB")


def main():
    if len(sys.argv) < 4:
        print("Usage: python3 extract_transcriptome.py <genome.fa> <annotation.gtf> <output.fa> [gene_name]")
        print("\nExample:")
        print("  python3 extract_transcriptome.py genome.fa annotation.gtf transcriptome.fa")
        print("  python3 extract_transcriptome.py genome.fa annotation.gtf CAMTA3_transcripts.fa CAMTA3")
        print("\nArguments:")
        print("  genome.fa      - Genome FASTA file (can be .gz)")
        print("  annotation.gtf - GTF annotation file (can be .gz)")
        print("  output.fa      - Output transcriptome FASTA file")
        print("  gene_name      - Optional: extract only this gene (for testing)")
        sys.exit(1)
    
    genome_file = sys.argv[1]
    gtf_file = sys.argv[2]
    output_file = sys.argv[3]
    gene_filter = sys.argv[4] if len(sys.argv) > 4 else None
    
    print("=" * 70)
    print("Transcriptome Extraction Tool (Fixed Version)")
    print("=" * 70)
    
    if gene_filter:
        print(f"\n⚠️  Gene filter enabled: Only extracting '{gene_filter}'")
    
    # 1. 解析 GTF
    try:
        transcripts = parse_gtf(gtf_file)
    except Exception as e:
        print(f"\n❌ ERROR parsing GTF file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    if len(transcripts) == 0:
        print("\n❌ ERROR: No transcripts found in GTF file!")
        print("Possible reasons:")
        print("  1. GTF file format is not standard")
        print("  2. No 'exon' or 'CDS' features in the file")
        print("  3. Missing transcript_id attribute")
        sys.exit(1)
    
    # 2. 加载基因组
    try:
        genome = load_genome(genome_file)
    except Exception as e:
        print(f"\n❌ ERROR loading genome file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    if len(genome) == 0:
        print("\n❌ ERROR: No sequences found in genome file!")
        sys.exit(1)
    
    # 3. 检查染色体名称匹配
    print("\nChecking chromosome name compatibility...")
    gtf_chrs = set(t['chr'] for t in transcripts.values() if t['chr'])
    genome_chrs = set(genome.keys())
    
    common_chrs = gtf_chrs & genome_chrs
    print(f"  GTF chromosomes: {len(gtf_chrs)}")
    print(f"  Genome chromosomes: {len(genome_chrs)}")
    print(f"  Matching chromosomes: {len(common_chrs)}")
    
    if len(common_chrs) == 0:
        print("\n⚠️  WARNING: No matching chromosome names!")
        print("  GTF chromosome examples:", list(gtf_chrs)[:5])
        print("  Genome chromosome examples:", list(genome_chrs)[:5])
        print("\n  This may cause extraction failures. Consider:")
        print("  1. Check if chromosome naming is consistent (e.g., '1' vs 'Chr1')")
        print("  2. Verify GTF and genome files are from the same assembly")
    
    # 4. 提取转录本序列
    try:
        write_transcriptome(transcripts, genome, output_file, gene_filter)
    except Exception as e:
        print(f"\n❌ ERROR extracting transcriptome: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print("\n" + "=" * 70)
    print("✓ Done!")
    print("=" * 70)


if __name__ == '__main__':
    main()
