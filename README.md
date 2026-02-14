# smRNA FISH Probe Design Pipeline

A complete pipeline for designing single molecule RNA FISH (smRNA FISH) probes with off-target checking.

## Overview

This pipeline designs RNA FISH probes by:
1. Scanning target gene sequence for potential probes using thermodynamic calculations
2. Filtering by GC content and PNAS rules
3. Removing repeats using RepeatMasker
4. Checking specificity against transcriptome using BLAST (optional)
5. Classifying probes into Exon and Junction types
6. Generating order lists for probe synthesis

## Table of Contents

1. [Probe Finding Principles](#1-probe-finding-principles)
2. [Filter Details](#2-filter-details)
3. [RepeatMasker Installation](#3-repeatmasker-installation)
4. [Specificity Check with BLAST](#4-specificity-check-with-blast)
5. [Probe Classification](#5-probe-classification)
6. [Usage](#usage)

---

## 1. Probe Finding Principles

### 1.1 Thermodynamic Calculation

Oligostan uses RNA thermodynamic properties to select optimal probe sequences. The melting temperature (Tm) is calculated using the nearest-neighbor model.

#### ΔG Calculation Formula

The free energy change (ΔG) at 37°C is calculated as:

```
ΔG = Σ(ΔG_dinucleotide) - (T × ΔS)
```

Where:
- ΔG_dinucleotide: thermodynamic values for each dinucleotide pair
- T: temperature (37°C = 310K)
- ΔS: entropy change

#### Nearest-Neighbor Parameters (SantaLucia, 1998)

| Dinucleotide | ΔG (kcal/mol) |
|-------------|---------------|
| AA/TT | -0.2 |
| AC/GT | -1.5 |
| AG/CT | -0.9 |
| AT/AT | -1.0 |
| CA/GT | -1.0 |
| CC/GG | -2.2 |
| CG/CG | -1.2 |
| CT/AG | -1.4 |
| GA/TC | -0.8 |
| GC/GC | -2.4 |
| GG/CC | -1.5 |
| GT/AC | -1.0 |
| TA/TA | -0.3 |
| TC/GA | -1.4 |
| TG/CA | -1.0 |
| TT/AA | -0.4 |

#### Probe Scoring

Probes are scored based on their ΔG values relative to the target ΔG (default: -33 kcal/mol):

```
Score = 1 - 0.1 × |ΔG_probe - ΔG_target|
```

Only probes with score ≥ 0.9 are retained.

#### Target ΔG Range

The pipeline scans for probes with ΔG values ranging from -36 to -28 kcal/mol (step: 0.5), ensuring probes have appropriate binding strength.

---

## 2. Filter Details

### 2.1 GC Content Filter

Probes must have GC content between specified range:

```r
MinGC <- 0.40    # 40%
MaxGC <- 0.60    # 60%
```

GC content calculation:
```
GC% = (G + C) / total_length
```

### 2.2 PNAS Rules

Based on the publication (Mueller et al., NAR 2013), five rules are applied:

| Rule | Description | Threshold |
|------|-------------|------------|
| 1 | A content | < 28% |
| 2 | No "AAAA" runs | - |
| 3 | C content | 22-28% |
| 4 | No "CCCC" runs in first 12nt | - |
| 5 | No ≥4 C in any 6nt window (first 12nt) | - |

Default: `PNASfilterOption <- c(1, 2, 4, 5)` (rules 3 is often too stringent)

### 2.3 RepeatMasker Filter

Probes with >10% repeat content are filtered out:

```r
MaskedFilter <- TRUE
MaxRepeatPercent <- 10%
```

### 2.4 Overlap Filter

Non-overlapping probes are selected with minimum spacing:

```r
DistanceMinInterSonde <- 2  # minimum 2bp spacing
```

---

## 3. RepeatMasker Installation

### 3.1 Download and Install

```bash
# 1. Download RepeatMasker
cd /usr/local
wget http://www.repeatmasker.org/RepeatMasker-open-4-X-XX.tar.gz
tar xzf RepeatMasker-open-4-X-XX.tar.gz

# 2. Download RMBlast (search engine)
wget http://www.repeatmasker.org/rmblast-2-X-XX.tar.gz
tar xzf rmblast-2-X-XX.tar.gz
cp -R rmblast-2-X-XX/* RepeatMasker/

# 3. Download TRF
wget http://tandem.bu.edu/trf/trf409.linux64
mv trf409.linux64 RepeatMasker/trf
chmod +x RepeatMasker/trf

# 4. Configure
cd RepeatMasker
perl ./configure
```

### 3.2 Species Selection

For Arabidopsis thaliana (Rosids clade):

```bash
# Using species name
RepeatMasker -species arabidopsis your_sequence.fa

# Using clade (for Rosids)
RepeatMasker -species rosids your_sequence.fa
```

Common species/groups:
- `arabidopsis` - Arabidopsis thaliana
- `rosids` - Rosids clade (includes Arabidopsis, many plants)
- `eudicots` - Eudicots
- `mammal` - Mammals
- `human` - Human

### 3.3 Disable RepeatMasker (Optional)

If RepeatMasker is not available, disable it in the R script:

```r
MaskedFilter <- FALSE
```

---

## 4. Specificity Check with BLAST

### 4.1 Why BLAST?

Even after filtering, probes may have off-target binding. BLAST ensures probe specificity against the transcriptome.

### 4.2 Install BLAST

```bash
# Using conda (recommended)
conda install -c bioconda blast

# Or download from NCBI
# https://blast.ncbi.nlm.nih.gov/doc/blast/download/
```

### 4.3 BLAST Parameters

The pipeline uses BLASTN with short sequence settings:

```bash
blastn -query probes.fasta \
       -subject transcriptome.fa \
       -task blastn-short \
       -outfmt '6 qseqid sseqid pident length qlen qcovs evalue' \
       -evalue 100 \
       -word_size 6 \
       -dust no \
       -soft_masking false
```

### 4.4 Specificity Criteria

| Parameter | Threshold |
|----------|-----------|
| Min Target Identity | ≥90% |
| Min Target Coverage | ≥85% |
| Max Off-target Score | <80% |

Off-target score = Identity × Coverage

---

## 5. Probe Classification

### 5.1 Exon vs Junction Probes

Probes are classified based on their genomic location:

- **Exon probes**: Located entirely within a single exon
- **Junction probes**: Span exon-exon boundaries

### 5.2 Detection Method

1. Extract mRNA reference sequence from genomic DNA
2. Identify exon junction positions
3. For each probe:
   - Map back to mRNA sequence
   - Check if probe spans a junction point

### 5.3 FLAP Sequences

Different FLAP sequences are used for multiplexed imaging:

| Type | FLAP Sequence | Use Case |
|------|---------------|-----------|
| X | `CCTCCTAAGTTTCGAGCTGGACTCAGTG` | Primary imaging |
| Y | `TTACACTCGGACCTCGTCGACATGCATT` | Secondary imaging |

---

## 6. Usage

### Two Workflow Modes

#### Mode 1: Simple (No Specificity Check)

```bash
# Edit configuration in smRNA_probe_design.r
# - FileNameFasta: your gene name
# - PathNameFasta: working directory

Rscript scripts/smRNA_probe_design.r
```

Results in `Probes_{GeneName}/*.csv`

#### Mode 2: With Specificity Check

```bash
# Step 1: Extract target gene transcriptome
python3 scripts/extract_transcriptome.py \
    /path/to/genome.fa \
    /path/to/annotation.gtf \
    target_transcriptome.fa \
    TARGET_GENE

# Step 2: Run probe design with specificity check
Rscript scripts/probe_design_with_transcriptome.r
```

### Configuration Parameters

Edit these in the R scripts:

```r
# Gene name
FileNameFasta <- 'CAMTA3'
FileNameGDNA  <- 'gCAMTA3.fa'

# Paths
PathNameFasta <- '/path/to/working/directory/'

# Probe parameters
TailleSondeMin <- 15
TailleSondeMax <- 31
MinGC <- 0.40
MaxGC <- 0.60
DistanceMinInterSonde <- 2

# RepeatMasker
MaskedFilter <- TRUE

# Specificity check
SpecificityCheck <- TRUE
TargetGenePattern <- 'CAMTA3|AT2G22300'
MinTargetIdentity <- 90
MinTargetCoverage <- 85
MaxOffTargetScore <- 80
```

## Output Files

| File | Description |
|------|-------------|
| `Probes_{Gene}/Order_List_GroupA_Exon_X.csv` | Exon probes (FLAP-X) |
| `Probes_{Gene}/Order_List_GroupB_Junc_Y.csv` | Junction probes (FLAP-Y) |
| `Probes_{Gene}/Specificity_Full_Report.csv` | Detailed specificity report |

## Requirements

- R 3.5+
- Python 3.6+
- RepeatMasker
- BLAST+ (for specificity check)
- R packages: seqinr, zoo

## Installation

```bash
# Install R packages
Rscript -e "install.packages(c('seqinr', 'zoo'), repos='http://cran.us.r-project.org')"

# Install RepeatMasker
# See Section 3 above

# Install BLAST
conda install -c bioconda blast
```

## Citation

If you use this pipeline, please cite:
- Mueller F, et al. (2013) smiFISH and FISH-quant. Nucleic Acids Research
- RepeatMasker: Smit AFA, Hubley R, Green P
- BLAST: Altschul SF, et al.

## License

MIT License

## Author

Zhihe Cai

## Contact

zh.cai@pku.edu.cn
