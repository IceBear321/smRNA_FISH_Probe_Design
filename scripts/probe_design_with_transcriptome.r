################################################################################
#      Probe Design with Automatic Transcriptome Extraction & Specificity     #
#      (Integrated solution for genome FASTA + GTF annotation)                 #
################################################################################

# ================= USER CONFIGURATION =================

# 1. 目标基因 (替换为你的基因名)
FileNameFasta <- 'test_gene' 
FileNameGDNA  <- 'genomic_test_gene.fa'
PathNameFasta <- '/data2/czh/tools/smRNA_probe_design/'

# 2. 探针参数
TailleSondeMin <- 15
TailleSondeMax <- 31
MinGC <- 0.40
MaxGC <- 0.60
DistanceMinInterSonde <- 2
PNASfilterOption <- c(1, 2, 4, 5)

# 3. RepeatMasker
MaskedFilter <- TRUE
RepeatMaskerCommand <- '/usr/local/RepeatMasker/RepeatMasker -species arabidopsis '

# 4. ⭐ 参考基因组和注释文件
GenomeFasta <- 'Arabidopsis_thaliana.TAIR10.dna.toplevel.transcripts.fa'
AnnotationGTF <- '/data/czh/reference_genome/Araport11_nonChr_nonTranspos_nonMt_nonPt.gtf'

# 5. ⭐ 特异性检查配置
SpecificityCheck <- TRUE

# 目标基因的多种可能命名（用于匹配转录本）
# 建议包含: 基因名、基因ID、转录本ID前缀
TargetGenePattern <- 'test_gene|ATXGXXXXX'

# 特异性阈值
MinTargetIdentity <- 90
MinTargetCoverage <- 85
MaxOffTargetScore <- 80

# 6. FLAP 序列
FLAP_X_SEQ <- "CCTCCTAAGTTTCGAGCTGGACTCAGTG" 
FLAP_Y_SEQ <- "TTACACTCGGACCTCGTCGACATGCATT" 

################################################################################
#                       CORE LOGIC                                             #
################################################################################

if (!require("seqinr")) install.packages("seqinr", repos="http://cran.us.r-project.org")
if (!require("zoo")) install.packages("zoo", repos="http://cran.us.r-project.org")
library(seqinr)
library(zoo)

setwd(PathNameFasta)
name_base <- strsplit(FileNameFasta, "[.]")[[1]][1]
FileNameOutput <- paste('Probes_', name_base, sep="")
if (!file.exists(FileNameOutput)){ dir.create(FileNameOutput) }

# --- 辅助函数 ---
ConvertRNASeq2DeltaGat37 <- function(RNASeq){
  NbBase = nchar(RNASeq)
  RNASeq <- toupper(RNASeq)
  DimVal  <- c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
  dGVal37 <- c(-0.2,-1.5,-0.9,-1.0,-1.0,-2.2,-1.2,-1.4,-0.8,-2.4,-1.5,-1.0,-0.3,-1.4,-1.0,-0.4)
  RNADimThermoTable <- data.frame(DimVal, dGVal37, stringsAsFactors=F)
  dim_vec <- substring(RNASeq, 1:(NbBase-1), 2:NbBase)
  dG <- RNADimThermoTable$dGVal37[match(dim_vec, RNADimThermoTable$DimVal)]
  return(data.frame(dim=dim_vec, dG=dG))
}

dGCalc.RNA.37 <- function(RNASeq, ProbeLength=31, SaltConc=0.115){
  if (nchar(RNASeq) < ProbeLength) return(NULL)
  RNASeqConv <- ConvertRNASeq2DeltaGat37(RNASeq=RNASeq)
  AlldG <- rollapply(RNASeqConv$dG, ProbeLength-1, sum)
  AlldG <- AlldG - ((log(SaltConc)*-0.175)-.2)
  return(AlldG)
}

dG37ScoreCalc <- function(ThedG37, DesireddG = -33){
  valtemp <- abs(ThedG37 - DesireddG)
  valnorm <- (-0.1 * valtemp) + 1
  return(valnorm)
}

isOk4GCFilter <- function(seq, min, max){
  comp <- count(s2c(tolower(seq)), 1)
  gc <- (comp['g'] + comp['c']) / sum(comp)
  return(gc >= min & gc <= max)
}

isOk4PNAS <- function(s, opts) {
  a <- (count(s2c(tolower(s)),1)['a']/nchar(s)) < 0.28
  aaaa <- !grepl("aaaa", tolower(s))
  c_count <- (count(s2c(tolower(s)),1)['c']/nchar(s))
  c <- c_count > 0.22 & c_count < 0.28
  cccc <- !grepl("cccc", tolower(s))
  sub <- substr(tolower(s), 1, 12)
  if(nchar(sub)<6) spec <- TRUE else spec <- !any(rollapply(s2c(sub), 6, function(x) sum(x=='c')) >= 4)
  res <- c(a, aaaa, c, cccc, spec)
  all(res[opts])
}

# ⭐ 新增：提取转录组序列
ExtractTranscriptome <- function(GenomeFasta, AnnotationGTF, TargetGene, OutputFasta) {
  cat("\n=== Extracting Transcriptome Sequences ===\n")
  cat(paste("Genome:", GenomeFasta, "\n"))
  cat(paste("Annotation:", AnnotationGTF, "\n"))
  cat(paste("Target gene:", TargetGene, "\n\n"))
  
  # 检查文件存在性
  if(!file.exists(GenomeFasta)) {
    cat("ERROR: Genome FASTA not found!\n")
    return(FALSE)
  }
  if(!file.exists(AnnotationGTF)) {
    cat("ERROR: GTF annotation not found!\n")
    return(FALSE)
  }
  
  # 查找 Python 提取脚本
  ScriptPath <- paste0(PathNameFasta, "extract_transcriptome.py")
  
  if(!file.exists(ScriptPath)) {
    cat("ERROR: extract_transcriptome.py not found!\n")
    cat(paste("Expected location:", ScriptPath, "\n"))
    cat("Please ensure the Python script is in the working directory.\n")
    return(FALSE)
  }
  
  # 运行提取脚本（只提取目标基因，加快速度）
  ExtractCmd <- paste0(
    "python3 ", ScriptPath, " ",
    GenomeFasta, " ",
    AnnotationGTF, " ",
    OutputFasta, " ",
    TargetGene
  )
  
  cat("Running transcriptome extraction...\n")
  ExitCode <- system(ExtractCmd, ignore.stdout = FALSE)
  
  if(ExitCode != 0) {
    cat("\nERROR: Transcriptome extraction failed!\n")
    return(FALSE)
  }
  
  if(!file.exists(OutputFasta) || file.info(OutputFasta)$size == 0) {
    cat("\nERROR: Output transcriptome file is empty or not created!\n")
    return(FALSE)
  }
  
  cat(paste("\n✓ Transcriptome extracted successfully:", OutputFasta, "\n"))
  return(TRUE)
}

# ⭐ 特异性检查函数
CheckProbeSpecificity <- function(ProbeDF, TargetGenePattern, TranscriptomeDB, 
                                   MinTargetIdent, MinTargetCov, MaxOffTargetScore) {
  
  cat("\n=== BLAST Specificity Check ===\n")
  
  # 创建临时 FASTA
  TmpProbeFile <- paste0(FileNameOutput, "/TempProbes_ForBLAST.fasta")
  write.fasta(as.list(ProbeDF$Seq), 
              names=as.character(1:nrow(ProbeDF)), 
              file.out=TmpProbeFile)
  
  # 运行 BLAST
  BlastOutFile <- paste0(FileNameOutput, "/blast_specificity_results.txt")
  BlastCmd <- paste0(
    "blastn",
    " -query ", TmpProbeFile,
    " -subject ", TranscriptomeDB,
    " -task blastn-short",
    " -outfmt '6 qseqid sseqid pident length qlen qcovs evalue'",
    " -evalue 100",
    " -word_size 6",
    " -dust no",
    " -soft_masking false",
    " -out ", BlastOutFile
  )
  
  cat("Running BLAST...\n")
  system(BlastCmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  # 检查输出
  if(!file.exists(BlastOutFile) || file.info(BlastOutFile)$size == 0) {
    cat("WARNING: BLAST returned no results\n")
    return(data.frame(ProbeDF, IsSpecific=FALSE, SpecDetails="No BLAST hits"))
  }
  
  # 解析结果
  BlastRes <- read.table(BlastOutFile, stringsAsFactors=F, sep="\t")
  colnames(BlastRes) <- c("ProbeID", "Target", "Identity", "Length", "QueryLen", "Coverage", "Evalue")
  
  cat(paste("Total BLAST hits:", nrow(BlastRes), "\n"))
  
  # 评估特异性
  IsSpecific <- logical(nrow(ProbeDF))
  SpecDetails <- character(nrow(ProbeDF))
  BestTarget <- character(nrow(ProbeDF))
  TargetIdentity <- numeric(nrow(ProbeDF))
  TargetCoverage <- numeric(nrow(ProbeDF))
  
  for(i in 1:nrow(ProbeDF)) {
    ProbeHits <- BlastRes[BlastRes$ProbeID == as.character(i), ]
    
    if(nrow(ProbeHits) == 0) {
      IsSpecific[i] <- FALSE
      SpecDetails[i] <- "No hits"
      next
    }
    
    ProbeHits$Score <- ProbeHits$Identity * ProbeHits$Coverage / 100
    ProbeHits <- ProbeHits[order(-ProbeHits$Score), ]
    
    # 查找目标基因
    TargetHits <- ProbeHits[grepl(TargetGenePattern, ProbeHits$Target, ignore.case=TRUE), ]
    
    HasGoodTarget <- FALSE
    if(nrow(TargetHits) > 0) {
      Best <- TargetHits[1, ]
      BestTarget[i] <- Best$Target
      TargetIdentity[i] <- Best$Identity
      TargetCoverage[i] <- Best$Coverage
      HasGoodTarget <- (Best$Identity >= MinTargetIdent) & (Best$Coverage >= MinTargetCov)
    } else {
      BestTarget[i] <- ProbeHits[1, "Target"]
      TargetIdentity[i] <- ProbeHits[1, "Identity"]
      TargetCoverage[i] <- ProbeHits[1, "Coverage"]
    }
    
    # 查找脱靶
    OffTargetHits <- ProbeHits[!grepl(TargetGenePattern, ProbeHits$Target, ignore.case=TRUE), ]
    HasOffTarget <- FALSE
    if(nrow(OffTargetHits) > 0 && OffTargetHits[1, "Score"] >= MaxOffTargetScore) {
      HasOffTarget <- TRUE
      SpecDetails[i] <- paste0("Off-target: ", OffTargetHits[1, "Target"])
    }
    
    IsSpecific[i] <- HasGoodTarget & !HasOffTarget
    if(IsSpecific[i]) SpecDetails[i] <- "Specific"
    else if(!HasGoodTarget) SpecDetails[i] <- "Poor target match"
  }
  
  ProbeDF$IsSpecific <- IsSpecific
  ProbeDF$SpecDetails <- SpecDetails
  ProbeDF$BestTarget <- BestTarget
  ProbeDF$TargetIdentity <- TargetIdentity
  ProbeDF$TargetCoverage <- TargetCoverage
  
  cat(paste("Specific probes:", sum(IsSpecific), "/", nrow(ProbeDF), "\n"))
  
  return(ProbeDF)
}

################################################################################
#                       MAIN WORKFLOW                                          #
################################################################################

# --- 0. 准备转录组数据库 ---
if(SpecificityCheck) {
  TranscriptomeDB <- paste0(FileNameOutput, "/", FileNameFasta, "_transcriptome.fa")
  
  if(!file.exists(TranscriptomeDB)) {
    cat("\n>>> [0/7] Preparing Transcriptome Database...\n")
    
    # 从基因名中提取主要部分（去掉可能的后缀）
    TargetGeneForExtraction <- strsplit(FileNameFasta, "[._]")[[1]][1]
    
    Success <- ExtractTranscriptome(
      GenomeFasta, 
      AnnotationGTF, 
      TargetGeneForExtraction, 
      TranscriptomeDB
    )
    
    if(!Success) {
      cat("\n❌ Failed to prepare transcriptome database!\n")
      cat("Options:\n")
      cat("  1. Fix the paths to genome and GTF files\n")
      cat("  2. Manually create transcriptome FASTA\n")
      cat("  3. Disable specificity check: SpecificityCheck <- FALSE\n\n")
      quit(save="no", status=1)
    }
  } else {
    cat(paste("\n>>> [0/7] Using existing transcriptome:", TranscriptomeDB, "\n"))
  }
}

# --- 1. 生成原始探针 ---
cat("\n>>> [1/7] Scanning sequence for probes...\n")
FullFastaPath <- paste0(PathNameFasta, FileNameFasta, ".fasta")
multifasta <- read.fasta(FullFastaPath)
TargetSeqBase <- paste(toupper(multifasta[[1]]), collapse="")

dG_Range <- seq(-36, -28, 0.5)
ProbesList <- list()

for(dG in dG_Range) {
  RawDG <- dGCalc.RNA.37(TargetSeqBase, ProbeLength=TailleSondeMax)
  if(is.null(RawDG)) next
  Scores <- dG37ScoreCalc(RawDG, dG)
  ValidIdx <- which(Scores >= 0.9)
  
  if(length(ValidIdx) > 0) {
    for(idx in ValidIdx) {
      StartPos <- idx
      EndPos <- idx + TailleSondeMax - 1
      TargetFragment <- substr(TargetSeqBase, StartPos, EndPos)
      ProbeSeq <- c2s(rev(comp(s2c(TargetFragment))))
      
      ProbesList[[length(ProbesList)+1]] <- data.frame(
        dGOpt = dG,
        ProbeName = paste0(FileNameFasta, "_dG", dG, "_", length(ProbesList)+1),
        StartPos = StartPos,
        Seq = ProbeSeq,
        dGVal = RawDG[idx],
        stringsAsFactors = F
      )
    }
  }
}

AllProbesDF <- do.call(rbind, ProbesList)
cat(paste("    Total Candidates:", nrow(AllProbesDF), "\n"))

AllProbesDF <- AllProbesDF[!duplicated(AllProbesDF$Seq), ]
cat(paste("    Unique Sequences:", nrow(AllProbesDF), "\n"))

# --- 2-3. 过滤 ---
cat("\n>>> [2/7] Applying GC & PNAS Filters...\n")
FilterResults <- logical(nrow(AllProbesDF))
for(i in 1:nrow(AllProbesDF)) {
  FilterResults[i] <- isOk4GCFilter(AllProbesDF$Seq[i], MinGC, MaxGC) & 
                      isOk4PNAS(AllProbesDF$Seq[i], PNASfilterOption)
}
FilteredDF <- AllProbesDF[FilterResults, ]
cat(paste("    Passing filters:", nrow(FilteredDF), "\n"))

if(MaskedFilter && nrow(FilteredDF) > 0) {
  cat("\n>>> [3/7] Running RepeatMasker...\n")
  TmpFa <- "RM_Input_Temp.fasta"
  write.fasta(as.list(FilteredDF$Seq), names=as.character(1:nrow(FilteredDF)), file.out=TmpFa)
  system(paste(RepeatMaskerCommand, TmpFa), ignore.stdout = TRUE)
  
  MaskedFile <- paste0(TmpFa, ".masked")
  if(file.exists(MaskedFile)) {
    MaskedSeqs <- read.fasta(MaskedFile)
    KeepRM <- logical(nrow(FilteredDF))
    for(i in 1:length(MaskedSeqs)) {
      s <- MaskedSeqs[[i]]
      n_count <- count(s, 1)['n']
      if(is.na(n_count)) n_count <- 0
      KeepRM[i] <- (n_count / length(s)) < 0.1
    }
    FilteredDF <- FilteredDF[KeepRM, ]
    cat(paste("    Passing RepeatMasker:", nrow(FilteredDF), "\n"))
    system(paste("rm -f", TmpFa, paste0(TmpFa, "*")), ignore.stderr = TRUE)
  }
} else {
  cat("\n>>> [3/7] Skipping RepeatMasker\n")
}

# --- 4. 特异性检查 ---
if(SpecificityCheck && nrow(FilteredDF) > 0) {
  cat("\n>>> [4/7] Checking Probe Specificity...\n")
  
  FilteredDF <- CheckProbeSpecificity(
    FilteredDF, 
    TargetGenePattern, 
    TranscriptomeDB,
    MinTargetIdentity,
    MinTargetCoverage,
    MaxOffTargetScore
  )
  
  # 保存完整报告
  FullReportFile <- paste0(FileNameOutput, "/Specificity_Full_Report.csv")
  write.csv(FilteredDF, FullReportFile, row.names=F, quote=F)
  cat(paste("Full report:", FullReportFile, "\n"))
  
  # 只保留特异性探针
  FilteredDF <- FilteredDF[FilteredDF$IsSpecific, ]
  
} else {
  cat("\n>>> [4/7] Skipping Specificity Check\n")
}

# --- 5-7. 去重叠、分类、输出 ---
cat("\n>>> [5/7] Removing Overlapping Probes...\n")
FinalProbes <- data.frame()

if(nrow(FilteredDF) > 0) {
  FilteredDF <- FilteredDF[order(FilteredDF$StartPos), ]
  LastEnd <- -999
  for(i in 1:nrow(FilteredDF)) {
    if(FilteredDF$StartPos[i] > (LastEnd + DistanceMinInterSonde)) {
      FinalProbes <- rbind(FinalProbes, FilteredDF[i, ])
      LastEnd <- FilteredDF$StartPos[i] + 32 - 1
    }
  }
}
cat(paste("    Final probes:", nrow(FinalProbes), "\n"))

if(nrow(FinalProbes) == 0) {
  cat("\n❌ No probes passed all filters!\n")
  quit(save="no", status=1)
}

cat("\n>>> [6/7] Classifying based on gDNA...\n")
gDNA_Path <- paste0(PathNameFasta, FileNameGDNA)
ProbeTypes <- rep("Exon (Unverified)", nrow(FinalProbes))

if(file.exists(gDNA_Path)) {
  gDNA_Lines <- readLines(gDNA_Path)
  gDNA_Seq <- paste(gDNA_Lines[!grepl("^>", gDNA_Lines)], collapse="")
  ExonMatches <- gregexpr("[A-Z]+", gDNA_Seq)[[1]]
  ExonLengths <- attr(ExonMatches, "match.length")
  
  mRNA_Ref_Seq <- ""
  Junction_Points <- c()
  Current_Len <- 0
  
  for(i in 1:length(ExonLengths)) {
    SeqFragment <- substr(gDNA_Seq, ExonMatches[i], ExonMatches[i] + ExonLengths[i] - 1)
    mRNA_Ref_Seq <- paste0(mRNA_Ref_Seq, SeqFragment)
    Current_Len <- Current_Len + nchar(SeqFragment)
    if(i < length(ExonLengths)) Junction_Points <- c(Junction_Points, Current_Len)
  }
  
  for(i in 1:nrow(FinalProbes)) {
    TargetFrag <- toupper(c2s(rev(comp(s2c(FinalProbes$Seq[i])))))
    MatchPos <- regexpr(TargetFrag, mRNA_Ref_Seq, fixed=TRUE)[1]
    
    if(MatchPos != -1) {
      MatchEnd <- MatchPos + nchar(TargetFrag) - 1
      IsJunction <- any(Junction_Points > MatchPos & Junction_Points < MatchEnd)
      ProbeTypes[i] <- ifelse(IsJunction, "Junction", "Exon")
    }
  }
}
FinalProbes$Type <- ProbeTypes
print(table(FinalProbes$Type))

cat("\n>>> [7/7] Generating Order CSVs...\n")
ExonList <- FinalProbes[grepl("Exon", FinalProbes$Type), ]
JuncList <- FinalProbes[FinalProbes$Type == "Junction", ]

if(nrow(ExonList) > 0) {
  OrderA <- data.frame(
    Name = paste0(FileNameFasta, "_Exon_X_", 1:nrow(ExonList)),
    Sequence = paste0(ExonList$Seq, FLAP_X_SEQ)
  )
  write.csv(OrderA, paste0(FileNameOutput, "/Order_List_GroupA_Exon_X.csv"), row.names=F, quote=F)
}

if(nrow(JuncList) > 0) {
  OrderB <- data.frame(
    Name = paste0(FileNameFasta, "_Junc_Y_", 1:nrow(JuncList)),
    Sequence = paste0(JuncList$Seq, FLAP_Y_SEQ)
  )
  write.csv(OrderB, paste0(FileNameOutput, "/Order_List_GroupB_Junc_Y.csv"), row.names=F, quote=F)
}

cat("\n")
cat("✓ COMPLETED\n")
cat(paste("Total probes:", nrow(FinalProbes), "\n"))
cat(paste("  - Exon:", nrow(ExonList), "\n"))
cat(paste("  - Junction:", nrow(JuncList), "\n"))
cat(paste("\nOutput folder:", FileNameOutput, "\n"))

