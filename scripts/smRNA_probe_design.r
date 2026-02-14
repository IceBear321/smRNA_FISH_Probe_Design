################################################################################
#                   Oligostan + Auto Classification + Order List               #
#                       (Fix: Loose Filter & Debug Log)                        #
################################################################################

# ================= USER CONFIGURATION (用户配置区) =================

# 1. 基因名
FileNameFasta <- 'CAMTA3' 
FileNameGDNA  <- 'gCAMTA3.fa'

# 2. 文件路径 (必须以 "/" 结尾)
PathNameFasta <- '/data2/czh/tools/smRNA_probe_design/'

# 3. 探针设计参数
TailleSondeMin <- 15
TailleSondeMax <- 31
MinGC <- 0.40
MaxGC <- 0.60
DistanceMinInterSonde <- 2 # 探针间最小间距
PNASfilterOption <- c(1, 2, 4, 5) 

# 4. RepeatMasker 设置
MaskedFilter <- TRUE
RepeatMaskerCommand <- '/usr/local/RepeatMasker/RepeatMasker -species arabidopsis ' 

# 5. FLAP 序列
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

ResultFileName <- paste(FileNameOutput, "/", FileNameOutput, "_FILT.txt", sep="")

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
  # Simplified PNAS check
  a <- (count(s2c(tolower(s)),1)['a']/nchar(s)) < 0.28
  aaaa <- !grepl("aaaa", tolower(s))
  c_count <- (count(s2c(tolower(s)),1)['c']/nchar(s))
  c <- c_count > 0.22 & c_count < 0.28
  cccc <- !grepl("cccc", tolower(s))
  # SpecStack simplified
  sub <- substr(tolower(s), 1, 12)
  if(nchar(sub)<6) spec <- TRUE else spec <- !any(rollapply(s2c(sub), 6, function(x) sum(x=='c')) >= 4)
  
  res <- c(a, aaaa, c, cccc, spec)
  all(res[opts])
}

# --- 1. 生成原始探针 ---
cat(">>> [1/6] Scanning sequence for probes...\n")
FullFastaPath <- paste0(PathNameFasta, FileNameFasta, ".fasta")
multifasta <- read.fasta(FullFastaPath)
TargetSeqBase <- paste(toupper(multifasta[[1]]), collapse="")

dG_Range <- seq(-36, -28, 0.5)
ProbesList <- list()

for(dG in dG_Range) {
  RawDG <- dGCalc.RNA.37(TargetSeqBase, ProbeLength=TailleSondeMax)
  if(is.null(RawDG)) next
  Scores <- dG37ScoreCalc(RawDG, dG)
  ValidIdx <- which(Scores >= 0.9) # Score cutoff
  
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
cat(paste("    Total Candidates Generated (Redundant):", nrow(AllProbesDF), "\n"))

# --- 2. 基础去重 (Same Seq) ---
# 保留 dG 分数最好的那个副本，或者是第一个
AllProbesDF <- AllProbesDF[!duplicated(AllProbesDF$Seq), ]
cat(paste("    Unique Sequences (After deduplication):", nrow(AllProbesDF), "\n"))

# --- 3. 过滤 (GC & PNAS) ---
cat(">>> [2/6] Applying GC & PNAS Filters...\n")
FilterResults <- logical(nrow(AllProbesDF))
for(i in 1:nrow(AllProbesDF)) {
  FilterResults[i] <- isOk4GCFilter(AllProbesDF$Seq[i], MinGC, MaxGC) & isOk4PNAS(AllProbesDF$Seq[i], PNASfilterOption)
}
FilteredDF <- AllProbesDF[FilterResults, ]
cat(paste("    Probes passing Filters:", nrow(FilteredDF), "\n"))

# --- 4. RepeatMasker ---
if(MaskedFilter && nrow(FilteredDF) > 0) {
  cat(">>> [3/6] Running RepeatMasker...\n")
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
    cat(paste("    Probes passing RepeatMasker:", nrow(FilteredDF), "\n"))
    system(paste("rm", TmpFa, paste0(TmpFa, "*")), ignore.stderr = T)
  }
}

# --- 5. 去重叠 (Greedy Selection) ---
cat(">>> [4/6] Removing Overlapping Probes (Spacing >", DistanceMinInterSonde, "nt)...\n")
FinalProbes <- data.frame()
if(nrow(FilteredDF) > 0) {
  FilteredDF <- FilteredDF[order(FilteredDF$StartPos), ]
  
  LastEnd <- -999
  CountDropped <- 0
  
  for(i in 1:nrow(FilteredDF)) {
    currStart <- FilteredDF$StartPos[i]
    if(currStart > (LastEnd + DistanceMinInterSonde)) {
      FinalProbes <- rbind(FinalProbes, FilteredDF[i, ])
      LastEnd <- FilteredDF$StartPos[i] + 32 - 1 # Approx length
    } else {
      CountDropped <- CountDropped + 1
    }
  }
  cat(paste("    Dropped due to overlap:", CountDropped, "\n"))
  cat(paste("    FINAL Non-overlapping Probes:", nrow(FinalProbes), "\n"))
}

# --- 6. gDNA 分类 ---
cat(">>> [5/6] Classifying based on gDNA...\n")
gDNA_Path <- paste0(PathNameFasta, FileNameGDNA)
ProbeTypes <- rep("Exon (Unverified)", nrow(FinalProbes)) # 默认值

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
table(FinalProbes$Type) # 打印分类结果

# --- 7. 生成订单 ---
cat(">>> [6/6] Generating Order CSVs...\n")
# 修正：包含 "Exon" 和 "Exon (Unverified)"
ExonList <- FinalProbes[grepl("Exon", FinalProbes$Type), ]
JuncList <- FinalProbes[FinalProbes$Type == "Junction", ]

# 输出 A
if(nrow(ExonList) > 0) {
  OrderA <- data.frame(Name = paste0(FileNameFasta, "_Exon_X_", 1:nrow(ExonList)),
                       Sequence = paste0(ExonList$Seq, FLAP_X_SEQ))
  write.csv(OrderA, paste0(FileNameOutput, "/Order_List_GroupA_Exon_X.csv"), row.names=F, quote=F)
}
# 输出 B
if(nrow(JuncList) > 0) {
  OrderB <- data.frame(Name = paste0(FileNameFasta, "_Junc_Y_", 1:nrow(JuncList)),
                       Sequence = paste0(JuncList$Seq, FLAP_Y_SEQ))
  write.csv(OrderB, paste0(FileNameOutput, "/Order_List_GroupB_Junc_Y.csv"), row.names=F, quote=F)
}
cat("Done. Check", FileNameOutput, "folder.\n")
