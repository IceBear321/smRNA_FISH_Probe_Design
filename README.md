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

1. [Environment Setup](#1-environment-setup)
2. [Probe Finding Principles](#2-probe-finding-principles)
3. [Filter Details](#3-filter-details)
4. [RepeatMasker Installation](#4-repeatmasker-installation)
5. [Specificity Check with BLAST](#5-specificity-check-with-blast)
6. [Probe Classification](#6-probe-classification)
7. [Usage](#7-usage)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Environment Setup

### 1.1 Required Software and Versions

This pipeline has been tested on the following environment:

| Software | Version | Installation Path |
|----------|---------|-------------------|
| R | 4.5.1 | `/usr/bin/R` |
| Python | 3.8.19 | `/home/czh/miniconda3/bin/python3` |
| BLAST+ | 2.12.0 | `/usr/bin/blastn` |
| RepeatMasker | 4.1+ | `/usr/local/RepeatMasker` |
| TRF | 4.09 | `/usr/local/bin/trf` |

### 1.2 R Packages Installation

Install required R packages:

```bash
# Start R
R

# Install packages (choose CRAN mirror when prompted)
install.packages("seqinr")
install.packages("zoo")

# Verify installation
library(seqinr)
library(zoo)
```

Or use command line:
```bash
Rscript -e "install.packages(c('seqinr', 'zoo'), repos='http://cran.us.r-project.org')"
```

### 1.3 Python Packages

The pipeline uses Python 3 with standard libraries. No additional packages required for basic functionality.

### 1.4 Working Directory Setup

Create a working directory for your probe design:

```bash
mkdir -p /path/to/your/project
cd /path/to/your/project

# Required input files:
# - test_gene.fasta         # Target gene DNA sequence (coding region only, no UTR)
# - genomic_test_gene.fa    # Genomic DNA with introns for classification
```

### 1.5 Input File Format

**test_gene.fasta** (Target gene, coding sequence only):
```
>test_gene
ATGCGTACGTTAGCTAGCTAGCTAG...
```

**genomic_test_gene.fa** (Genomic DNA with introns):
```
>genomic_test_gene
ATGCGT... intron1 ... GTAG... intron2 ... ATGC...
```

**Important Notes:**
- Use only the CODING SEQUENCE (CDS) for probe design, exclude UTRs
- If CDS is very short (<500bp), UTRs can be included
- Do not use underscores in file names or sequence headers

---

## 2. Probe Finding Principles

### 2.1 Thermodynamic Calculation

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

#### Probe Length

- Minimum: 15-31 nucleotides (configurable)
- Default: 15-31 bp
- Longer probes may have better specificity but lower efficiency

---

## 3. Filter Details

### 3.1 GC Content Filter

Probes must have GC content between specified range:

```r
MinGC <- 0.40    # 40%
MaxGC <- 0.60    # 60%
```

GC content calculation:
```
GC% = (G + C) / total_length
```

**Recommendation:** GC content affects probe binding efficiency. 40-60% is optimal for smFISH.

### 3.2 PNAS Rules

Based on the publication (Mueller et al., NAR 2013), five rules are applied to ensure probe quality:

| Rule | Description | Threshold |
|------|-------------|------------|
| 1 | A (adenine) content | < 28% |
| 2 | No "AAAA" runs | - |
| 3 | C (cytosine) content | 22-28% |
| 4 | No "CCCC" runs in first 12nt | - |
| 5 | No ≥4 C in any 6nt window (first 12nt) | - |

Default filter option: `PNASfilterOption <- c(1, 2, 4, 5)`

**Note:** Rule 3 is often too stringent and may be excluded. If too few probes are obtained, try removing rules gradually.

### 3.3 RepeatMasker Filter

Probes with >10% repeat content are filtered out:

```r
MaskedFilter <- TRUE
MaxRepeatPercent <- 10%
```

RepeatMasker identifies and masks:
- Transposable elements
- Simple repeats (microsatellites)
- Low complexity regions

### 3.4 Overlap Filter (Greedy Selection)

Non-overlapping probes are selected with minimum spacing:

```r
DistanceMinInterSonde <- 2  # minimum 2bp spacing between probes
```

Algorithm:
1. Sort probes by genomic position
2. Select first probe
3. Skip overlapping probes until minimum spacing is met
4. Repeat until all probes processed

---

## 4. RepeatMasker Installation

### 4.1 Download Components

Download from official websites:

```bash
# Create installation directory
sudo mkdir -p /usr/local/RepeatMasker
cd /usr/local

# 1. RepeatMasker core
wget http://www.repeatmasker.org/RepeatMasker-open-4-X-XX.tar.gz
tar xzf RepeatMasker-open-4-X-XX.tar.gz

# 2. RMBlast (required search engine)
wget http://www.repeatmasker.org/rmblast-2-X-XX.tar.gz
tar xzf rmblast-2-X-XX.tar.gz
cp -R rmblast-2-X-XX/* RepeatMasker/

# 3. TRF (Tandem Repeat Finder)
wget http://tandem.bu.edu/trf/trf409.linux64
mv trf409.linux64 RepeatMasker/trf
chmod +x RepeatMasker/trf

# 4. Repeat Library (Dfam or RepBase)
# For Arabidopsis (Rosids), recommend RepBase
# Register at https://girinst.org/ and download repeat library
# Or use Dfam for quick setup (limited species)
```

### 4.2 Configure RepeatMasker

```bash
cd /usr/local/RepeatMasker
perl ./configure
```

Configuration prompts:
```
Where is Perl? [/usr/bin/perl]
Where is RepeatMasker? [/usr/local/RepeatMasker]
Where is TRF? [/usr/local/RepeatMasker]
Search engine [2]: 2  # RMblast
RMblast directory [/usr/local/RepeatMasker/rmblast-X.X.X]
```

### 4.3 Species Selection for Arabidopsis

For Arabidopsis thaliana (or other Rosids):

```bash
# Using species name
RepeatMasker -species arabidopsis your_sequence.fa

# Using clade (broader coverage)
RepeatMasker -species rosids your_sequence.fa

# Common plant species:
# - arabidopsis
# - rosids (includes Arabidopsis, Brassica, etc.)
# - eudicots
# - viridiplantae (all green plants)
```

**In the R script, modify:**
```r
RepeatMaskerCommand <- '/usr/local/RepeatMasker/RepeatMasker -species arabidopsis '
```

### 4.4 Disable RepeatMasker (Optional)

If RepeatMasker is not available or not working:

```r
# In the R script, change:
MaskedFilter <- FALSE
```

---

## 5. Specificity Check with BLAST

### 5.1 Why BLAST?

Even after filtering, probes may have off-target binding. BLAST ensures probe specificity against the transcriptome.

### 5.2 Install BLAST

```bash
# Using conda (recommended)
conda install -c bioconda blast

# Or from source
# Download from: https://blast.ncbi.nlm.nih.gov/doc/blast/download/
```

Verify installation:
```bash
blastn -version
# Should show: blastn: 2.12.0+
```

### 5.3 Prepare Transcriptome Database

For specificity checking, you need a transcriptome FASTA file:

```bash
# Option 1: Download from Ensembl/NCBI
# Arabidopsis: https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index

# Option 2: Extract from genome using provided script
python3 scripts/extract_transcriptome.py \
    /path/to/genome.fa \
    /path/to/annotation.gtf \
    transcriptome.fa \
    test_gene
```

### 5.4 BLAST Parameters

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

### 5.5 Specificity Criteria

| Parameter | Threshold | Description |
|----------|-----------|-------------|
| Min Target Identity | ≥90% | Probe matches target gene |
| Min Target Coverage | ≥85% | Most of probe binds to target |
| Max Off-target Score | <80% | No strong off-target binding |

Off-target score calculation:
```
Score = Identity × Coverage
```

---

## 6. Probe Classification

### 6.1 Exon vs Junction Probes

Probes are classified based on their genomic location:

- **Exon probes**: Located entirely within a single exon
- **Junction probes**: Span exon-exon boundaries

### 6.2 Detection Method

1. Extract mRNA reference sequence from genomic DNA
2. Identify exon junction positions (splice sites)
3. For each probe:
   - Reverse complement (probe is RNA, genomic is DNA)
   - Map back to mRNA sequence
   - Check if probe spans a junction point (exists between exons)

### 6.3 FLAP Sequences

Different FLAP sequences are used for multiplexed imaging:

| Type | FLAP Sequence | Use Case |
|------|---------------|----------|
| X | `CCTCCTAAGTTTCGAGCTGGACTCAGTG` | Primary imaging channel |
| Y | `TTACTCGGACCTCGTCGACATGCATT` | Secondary imaging channel |

**Note:** These are concatenated to the 3' end of probe sequences for ordering.

---

## 7. Usage

### 7.1 Quick Start

```bash
# 1. Prepare working directory
mkdir my_probe_design
cd my_probe_design

# 2. Create input files
# - test_gene.fasta (your target gene)
# - genomic_test_gene.fa (genomic sequence with introns)

# 3. Run simple mode (no specificity check)
Rscript /path/to/smRNA_probe_design.r
```

### 7.2 Configuration

Edit the following in the R script before running:

```r
# ================= USER CONFIGURATION =================

# 1. Gene name (replace test_gene with your gene name)
FileNameFasta <- 'test_gene'
FileNameGDNA  <- 'genomic_test_gene.fa'

# 2. File paths (must end with "/")
PathNameFasta <- '/path/to/your/working/directory/'

# 3. Probe design parameters
TailleSondeMin <- 15        # Minimum probe length
TailleSondeMax <- 31        # Maximum probe length
MinGC <- 0.40              # Minimum GC content
MaxGC <- 0.60              # Maximum GC content
DistanceMinInterSonde <- 2  # Minimum spacing between probes

# 4. PNAS filter options (which rules to apply)
PNASfilterOption <- c(1, 2, 4, 5)

# 5. RepeatMasker settings
MaskedFilter <- TRUE
RepeatMaskerCommand <- '/usr/local/RepeatMasker/RepeatMasker -species arabidopsis '

# 6. FLAP sequences (do not modify unless necessary)
FLAP_X_SEQ <- "CCTCCTAAGTTTCGAGCTGGACTCAGTG"
FLAP_Y_SEQ <- "TTACACTCGGACCTCGTCGACATGCATT"
```

### 7.3 Two Workflow Modes

#### Mode 1: Simple (No Specificity Check)

```bash
# Run basic pipeline
Rscript scripts/smRNA_probe_design.r
```

Results in `Probes_test_gene/*.csv`:
- `Order_List_GroupA_Exon_X.csv` - Exon probes with FLAP-X
- `Order_List_GroupB_Junc_Y.csv` - Junction probes with FLAP-Y

#### Mode 2: With Specificity Check

```bash
# Step 1: Extract target gene transcriptome (if needed)
python3 scripts/extract_transcriptome.py \
    /path/to/genome.fa \
    /path/to/annotation.gtf \
    test_gene_transcriptome.fa \
    test_gene

# Step 2: Run probe design with specificity check
Rscript scripts/probe_design_with_transcriptome.r

# Additional configuration for Mode 2:
# TargetGenePattern <- 'test_gene|ATXGXXXXX'  # Gene name and ID
# MinTargetIdentity <- 90
# MinTargetCoverage <- 85
# MaxOffTargetScore <- 80
```

---

## 8. Troubleshooting

### 8.1 No Probes Found

**Problem:** Output file is empty

**Solutions:**
1. Reduce PNAS filter stringency:
   ```r
   PNASfilterOption <- c(1, 2, 4)  # Remove rule 5
   PNASfilterOption <- c(1, 2)       # Even fewer rules
   ```

2. Widen GC content range:
   ```r
   MinGC <- 0.30
   MaxGC <- 0.70
   ```

3. Extend probe length range:
   ```r
   TailleSondeMin <- 12
   TailleSondeMax <- 35
   ```

4. Disable RepeatMasker:
   ```r
   MaskedFilter <- FALSE
   ```

### 8.2 RepeatMasker Not Working

**Problem:** "RepeatMasker: command not found"

**Solutions:**
1. Verify installation:
   ```bash
   which RepeatMasker
   /usr/local/RepeatMasker/RepeatMasker
   ```

2. Check path in R script matches your installation

3. Try disabling if still not working:
   ```r
   MaskedFilter <- FALSE
   ```

### 8.3 BLAST Specificity Check Fails

**Problem:** "blastn: command not found"

**Solutions:**
1. Install BLAST:
   ```bash
   conda install -c bioconda blast
   ```

2. Or disable specificity check:
   ```r
   SpecificityCheck <- FALSE
   ```

### 8.4 Low Specificity

**Problem:** Many probes have off-target binding

**Solutions:**
1. Increase specificity thresholds:
   ```r
   MinTargetIdentity <- 95
   MinTargetCoverage <- 90
   MaxOffTargetScore <- 50
   ```

2. Review Specificity_Full_Report.csv and manually filter

---

## Output Files

| File | Description |
|------|-------------|
| `Probes_{Gene}/Order_List_GroupA_Exon_X.csv` | Exon probes (FLAP-X) for ordering |
| `Probes_{Gene}/Order_List_GroupB_Junc_Y.csv` | Junction probes (FLAP-Y) for ordering |
| `Probes_{Gene}/Specificity_Full_Report.csv` | Detailed specificity report (Mode 2) |
| `Probes_{Gene}/{Gene}_FILT.txt` | Filtered probes with all details |

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
