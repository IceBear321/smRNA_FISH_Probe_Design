# smRNA FISH Probe Design Pipeline

A complete pipeline for designing single molecule RNA FISH (smRNA FISH) probes with off-target checking.

## Overview

This pipeline designs RNA FISH probes by:
1. Scanning target gene sequence for potential probes
2. Filtering by GC content and thermodynamic properties
3. Removing repeats using RepeatMasker
4. Checking specificity against transcriptome (optional)
5. Classifying probes into Exon and Junction types
6. Generating order lists for probe synthesis

## Two Workflow Modes

### Mode 1: Simple (No Specificity Check)
```bash
Rscript scripts/smRNA_probe_design.r
```
- Faster, suitable for initial probe design
- Results in `Probes_CAMTA3/*.csv`

### Mode 2: With Specificity Check
```bash
# Step 1: Extract transcriptome (if not exists)
python3 scripts/extract_transcriptome.py \
    /path/to/genome.fa \
    /path/to/annotation.gtf \
    transcriptome.fa \
    TARGET_GENE

# Step 2: Run probe design with specificity check
Rscript scripts/probe_design_with_transcriptome.r
```

## Required Data Files

| File | Description |
|------|-------------|
| `{GeneName}.fasta` | Target gene DNA sequence |
| `g{GeneName}.fa` | Genomic DNA with introns |
| `Arabidopsis_thaliana.TAIR10.dna.toplevel.transcripts.fa` | Reference transcriptome (for specificity check) |

## Output Files

| File | Description |
|------|-------------|
| `Probes_{GeneName}/Order_List_GroupA_Exon_X.csv` | Exon probes (FLAP-X) |
| `Probes_{GeneName}/Order_List_GroupB_Junc_Y.csv` | Junction probes (FLAP-Y) |
| `Probes_{GeneName}/Specificity_Full_Report.csv` | Detailed specificity report (if enabled) |

## Configuration Parameters

In the R scripts, modify these variables:

```r
# 1. Gene name
FileNameFasta <- 'CAMTA3'
FileNameGDNA  <- 'gCAMTA3.fa'

# 2. Paths
PathNameFasta <- '/path/to/working/directory/'

# 3. Probe parameters
TailleSondeMin <- 15
TailleSondeMax <- 31
MinGC <- 0.40
MaxGC <- 0.60

# 4. RepeatMasker (optional)
MaskedFilter <- TRUE

# 5. FLAP sequences
FLAP_X_SEQ <- "CCTCCTAAGTTTCGAGCTGGACTCAGTG"
FLAP_Y_SEQ <- "TTACACTCGGACCTCGTCGACATGCATT"
```

## Requirements

- R 3.5+
- Python 3.6+
- RepeatMasker
- BLAST (for specificity check)
- R packages: seqinr, zoo

## Installation

```bash
# Install R packages
Rscript -e "install.packages(c('seqinr', 'zoo'), repos='http://cran.us.r-project.org')"

# Install RepeatMasker
# See: http://www.repeatmasker.org/

# Install BLAST
conda install -c bioconda blast
```

## Probe Classification

- **Exon probes**: Located within single exons
- **Junction probes**: Span exon-exon junctions

Different FLAP sequences (X or Y) are used for multiplexed imaging.

## Citation

If you use this pipeline, please cite:
- smRNA FISH methodology
- RepeatMasker
- BLAST

## License

MIT License

## Author

Zhihe Cai

## Contact

zh.cai@pku.edu.cn
