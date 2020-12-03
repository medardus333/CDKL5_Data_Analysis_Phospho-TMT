## # # # # # # # # # # # # # # # # # # # # # # # #
##
## Project: CDKL5-NLS Phosphoproteomics
##
## Script name: Script_File_S2_Data_analysis_Phospho-TMT.R
##
## Purpose of script: Process MaxQuant output and data analysis
##
## Author: Florian Weiland
##
## Date Created: 2020-04-03
##
## # # # # # # # # # # # # # # # # # # # # # # # #

## Libraries

library(ggplot2)
library(reshape2)
library(vsn)
library(limma)
library(seqinr)
library(plyr)
library(stringr)
library(ggrepel)
library(ggpointdensity)
library(wesanderson)
library(extrafont)
library(scales)

## Functions ----
source("./Functions/averageMaxQuant.R")

## Set fraction numbers and reporter intensity groups ----
## Match the way the data in the "Experiment" column of the MaxQuant "evidence.txt" file is labelled
fractions <- seq(1,20,1)

## WT and KD columns in MaxQuant "evidence.txt" file
wt.cols <- "^Reporter.intensity.corrected\\.[1,2,3,4,5]$"
kd.cols <- "^Reporter.intensity.corrected.[6,7,8,9]|Reporter.intensity.corrected.10$"

## Set folders ----

## Where is MaxQuant "evidence.txt" stored?
mq.dir <- "./Data/"

## Where are databases for metadata stored?
db.dir <- "./Databases/"

## Where do results go?
output.res.dir <- "./Output_R/"

### Set file names ----

## All TMT results
results.file <- "Phospho-TMT_results_complete.csv"

## Cleaned up TMT results
results.clean.file <- "Phospho-TMT_results_clean.csv"

### Plots and images

## VSN fit
vsn.fit.file <- "MeanSDplot_VSN_fit.png"

## QQ-plot
qqplot.file <- "QQ-plot.png"

## VSN calibration check
VSN.transform.TMT.file <- "VSN_calibration.png"

## Garden Sprinkler plot
gs.file <- "Garden_sprinkler_plot.png"

## Unlabelled Garden Sprinkler plot
gs.file.2 <- "Garden_sprinkler_plot_unlabelled.png"

## Descriptive statistics

stats.file <- "Descriptive_Statistics.csv"

## Set ggplot theme and palette ----

theme_set(theme_classic(base_family = "Arial"))
pal <- wes_palette("Darjeeling1")

# # # # # # # # # # # # # # # # # # # # # # # 

## Read in data ----
message(paste0("\n\nLoading ", mq.dir, "evidence.txt"))

data.mq <- read.csv(
  paste0(mq.dir, "evidence.txt"),
  sep = "\t",
  stringsAsFactors = FALSE
)


## Remove Reverse hits, potential contaminants, PIF < 0.65 and PEP >= 0.05

rem.rev <- which(data.mq$Reverse == "+")
rem.cont <- which(data.mq$Potential.contaminant == "+")
rem.pep <- which(data.mq$PEP >= 0.05)
rem.pif <- which(data.mq$PIF < 0.60)
rem.pif.2 <- which(is.nan(data.mq$PIF) == TRUE)

data.preclean <- data.mq[-c(rem.rev, rem.cont, rem.pep, rem.pif, rem.pif.2), ]

## Count zero intensities for each quantified peptide, remove KD or WT rows with less than 2 datapoints

group.wt <- grep(
  wt.cols,
  colnames(data.preclean)
)

group.kd <- grep(
  kd.cols,
  colnames(data.preclean)
)

zero.kd <- rowSums(data.preclean[, group.kd] == 0)
zero.wt <- rowSums(data.preclean[, group.wt] == 0)

## Remove rows with >= 4 counts of 0 (= at least 2 datapoints in KD/WT group)

rem.zero.kd <- which(zero.kd >= 4)
rem.zero.wt <- which(zero.wt >= 4)

data.clean <- data.preclean[-c(rem.zero.kd, rem.zero.wt), ]

## Replace zeroes with NA in reporter intensities columns (MaxQuant puts zeroes instead NAs)

data.clean[, c(group.wt, group.kd)][data.clean[, c(group.wt, group.kd)] == 0] <- NA

### Data averaging ----

## Peptides with one observatiton are kept
## Peptides with two observations get averaged
## Peptide with three or more observation, peptides are averaged or peptide with median ratio WT/KD is taken
## This is being done for each fraction independently to avoid bias by fraction dependent co-isolation populations 

message("Averaging data")

data.av <- averageMaxQuant(data.clean, fractions)

## Data analysis ----
## VSN tansformation ----

message("\nVSN Transform")

data.av.matrix <- as.matrix(data.av[ , c(group.wt, group.kd)])

## colMeans function used in averageMaxQuant returns NaN
## if all observations of a peptide in a TMT channel are NA.
## Replace with NA, VSN cannot handle otherwise
## Also is.nan() does not work in dataframes...

data.av.matrix[is.nan(data.av.matrix)] <- NA 

vsn.model <- vsn2(data.av.matrix)

## Check model fit

message(paste0("Writing ", output.res.dir, vsn.fit.file))

png(
  paste0(output.res.dir, vsn.fit.file),
  width = 6,
  height = 6,
  unit = "in",
  res = 600
)

vsn.plot <- meanSdPlot(vsn.model)

vsn.plot$gg +
  ggtitle("VSN fit of TMT data")

dev.off()

## VSN transform TMT data

data.trans <- predict(vsn.model, newdata = data.av.matrix)

## QQ plot to check normality

message(paste0("Writing ", output.res.dir, qqplot.file))

png(
  paste0(output.res.dir, qqplot.file),
  width = 6,
  height = 6,
  unit = "in",
  res = 600
)

par(mfrow = c(1,2))

qqnorm(
  data.av.matrix,
  main = "Normal Q-Q- Plot \n TMT raw data"
)
qqline(data.av.matrix)

qqnorm(
  data.trans,
  main = "Normal Q-Q- Plot \n TMT VSN transformed data"
)
qqline(data.trans)

dev.off()

## Boxplot VSN intensities to check intensity calibration

message(paste0("Writing ", output.res.dir, VSN.transform.TMT.file))

data.box <- melt(data.trans, id = colnames(data.trans))
colnames(data.box) <- c("Index", "TMT.label", "VSN.Intensity")

data.box$Group <- data.box$TMT.label
data.box$Group <- sub(wt.cols, "WT", data.box$Group)
data.box$Group <- sub(kd.cols, "KD", data.box$Group)
data.box$Group <- factor(data.box$Group, levels = c("WT", "KD"))

data.box$TMT.label <- sub("Reporter.intensity.corrected.", "", data.box$TMT.label)

tmt.labels <- c(
  "126",
  "127N",
  "127C",
  "128N",
  "128C",
  "129N",
  "129C",
  "130N",
  "130C",
  "131"
)

reporter <- seq(1,10,1)

for (cur.tmt in 1:length(tmt.labels)) {

data.box$TMT.label <- sub(paste0("^", reporter[cur.tmt], "$"), tmt.labels[cur.tmt], data.box$TMT.label)

}

data.box$TMT.label <- factor(data.box$TMT.label, levels = tmt.labels)

ggplot(data.box, aes(x = TMT.label, y = VSN.Intensity, fill = Group)) +
  geom_boxplot(
    lwd = 0.8,
    na.rm = TRUE
  ) +
  scale_fill_manual(values = c(pal[1], pal[5])) +
  theme(
    axis.text.x = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 16),
    axis.text.y = element_text(family = "Roboto Mono", size = 14),
    axis.line = element_line(
		  color = "black", 
      size = 1,
		  linetype = "solid"
	  ),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(.20, "cm")
  ) +
  labs(
    x = "TMT label",
    y = "VSN transformed Reporter.intensity.corrected"
  )

ggsave(
  paste0(output.res.dir, VSN.transform.TMT.file),
  dpi = 600
)

## Limma ----
## Build dataframe for limma
## KD as base, goes therefore first in data.frame

message("Limma calculation")

data.limma <- data.frame(
  Sequence = data.av$Sequence,
  Modified.sequence = data.av$Modified.sequence,
  Phospho.STY.probabilities = data.av$Phospho..STY..Probabilities,
  Phospho.STY = data.av$Phospho..STY.,
  Leading.razor.protein = data.av$Leading.razor.protein,
  Gene.names =  data.av$Gene.names,
  Fraction =  data.av$Experiment,
  data.trans[, c(6:10, 1:5)],
  stringsAsFactors = FALSE
)

## Remove non-phospho peptides

rem.non.phos <- which(data.limma[, "Phospho.STY"] == 0)

data.limma <- data.limma[-rem.non.phos, ]

## Implement counter to be sure limma output can be correctly merged with MaxQuant data

data.limma$Counter <- seq(1, nrow(data.limma), 1)
rownames(data.limma) <- data.limma$Counter

## Get Reporter Intensity columns

int.limma <- grep("Reporter", colnames(data.limma))

## Build model matrix

groups <- c(rep("KD",5), rep("WT",5)) 
groups <- factor(groups, levels = c("KD", "WT"))

design <- model.matrix(~ groups) 

## Limma
## eBayes with robust = TRUE to protect from hypervariable observations

model.limma <- lmFit(data.limma[, int.limma], design) 
fit <- eBayes(model.limma, robust = TRUE) 

TMT.results <- topTable(fit, number = nrow(data.limma))

## Combine limma dataframe with TMT data

TMT.results <- TMT.results[order(as.numeric(rownames(TMT.results))), ]
TMT.results[, "Counter"] <- rownames(TMT.results)

TMT.results <- merge(data.limma, TMT.results, by = "Counter")

## Add Metadata to each peptide ----
## Read in databases

message(paste0("Loading databases for metadata from ", db.dir))

phos.site.plus <- read.csv(
  paste0(db.dir, "Phosphorylation_site_dataset.txt"),
  skip = 3,
  sep = "\t",
  stringsAsFactors = FALSE
)

uniprot.go <- read.csv(
  paste0(db.dir, "200405_Uniprot_human_GO.tab"),
  sep = "\t",
  stringsAsFactors = FALSE
)

uniprot.phos <- read.fasta(
  file = paste0(db.dir, "180405_Uniprot_Homo+sapiens_canonical+isoforms.fasta"),
  seqtype = "AA",
  strip.desc = TRUE,
  as.string = TRUE
)

## Extract protein annotations , AA sequence and combine
annot <- lapply(uniprot.phos, attributes)
annot <- unlist(annot)
annot <- as.vector(annot[seq(2, length(annot), 3)])

uniprot.phos <- as.data.frame(ldply(uniprot.phos, rbind))
colnames(uniprot.phos) <- c("Protein", "AA.seq")
uniprot.phos$AA.seq <- as.character(uniprot.phos$AA.seq)

uniprot.phos <- cbind(uniprot.phos, annot)
uniprot.phos$annot <- as.character(uniprot.phos$annot)

## Extract Gene names
uniprot.phos$Gene.names <- str_extract(uniprot.phos$annot, "GN\\=.*? |GN\\=.*?$")
uniprot.phos$Gene.names <- gsub("GN=| ", "", uniprot.phos$Gene.names)

## Extract Uniprot Acc. Number

uniprot.phos$Uniprot.Acc <- str_extract(uniprot.phos$annot, "\\|.*?\\|")
uniprot.phos$Uniprot.Acc <- gsub("\\|", "", uniprot.phos$Uniprot.Acc)

## Extract Protein name

uniprot.phos$protein.name.long <- str_extract(uniprot.phos$annot, " .*? OS\\=")
uniprot.phos$protein.name.long <- gsub("^ | OS\\=", "", uniprot.phos$protein.name.long)

## Check if peptide is isoform specific, if not use canonical form for p-site number.
## Absolute position of phospho-peptide in AA seq
## MaxQuant decides for some reason that peptide is part of specific isoform of protein even when found in
## canonical form. We want to have position of p-site in canonical form.

TMT.results$Uniprot.Acc <- TMT.results$Leading.razor.protein
TMT.results$Uniprot.Acc <- gsub("-\\d{1,2}", "", TMT.results$Uniprot.Acc)

gene.cols <- which(colnames(uniprot.phos) == "Gene.names")
uni.merge <- uniprot.phos[, -gene.cols]

TMT.results <- merge(TMT.results, uni.merge, by = "Uniprot.Acc")

pep.pos <- str_locate(TMT.results$AA.seq, TMT.results$Sequence)
TMT.results[, "Pep.pos.in.AA.seq"] <- pep.pos[, "start"]

iso.pep <- which(is.na(TMT.results$Pep.pos.in.AA.seq) == TRUE)

TMT.results$Isoform.specific.peptide <- FALSE
TMT.results[iso.pep, "Isoform.specific.peptide"] <- TRUE

## Put isoform specific Uniprot Accession number, protein name and AA.seq where iso.pep indicates

TMT.results[iso.pep, "Uniprot.Acc"] <- TMT.results[iso.pep, "Leading.razor.protein"]

for (cur.iso.pep in iso.pep) {
  
  temp.row.uni <- which(TMT.results[cur.iso.pep, "Uniprot.Acc"] == uniprot.phos$Uniprot.Acc)
  TMT.results[cur.iso.pep, "AA.seq"] <- uniprot.phos[temp.row.uni, "AA.seq"]
  TMT.results[cur.iso.pep, "protein.name.long"] <- uniprot.phos[temp.row.uni, "protein.name.long"]
  
}

pep.pos <- str_locate(TMT.results$AA.seq, TMT.results$Sequence)
TMT.results[, "Pep.pos.in.AA.seq"] <- pep.pos[, "start"]

## Phospho-Site plus uses non-isoform entries, need to de-isoform our protein names
TMT.results$Parent.Entry <- TMT.results$Leading.razor.protein
TMT.results$Parent.Entry <- gsub("-\\d{1,2}", "", TMT.results$Parent.Entry)

## Extract phospho-peptide sequence without indication of other modifications
TMT.results[, "Modified.sequence.2"] <- gsub("\\(ox\\)", "", TMT.results[, "Modified.sequence"])
TMT.results[, "Modified.sequence.2"] <- gsub("\\(ac\\)", "", TMT.results[, "Modified.sequence.2"])
TMT.results[, "Modified.sequence.2"] <- gsub("\\(de\\)", "", TMT.results[, "Modified.sequence.2"])
TMT.results[, "Modified.sequence.2"] <- gsub("\\(ca\\)", "", TMT.results[, "Modified.sequence.2"])
TMT.results[, "Modified.sequence.2"] <- gsub("_", "", TMT.results[, "Modified.sequence.2"])

## Add rest of metadata: Position of p-site and is phosphorylation site part of RPX-S/T motif?

message("Adding metadata")

progress.bar.2 <- txtProgressBar(min = 0, max = nrow(TMT.results), style = 3)

pep.final <- data.frame()

TMT.meta <- TMT.results[, c("Counter", "Modified.sequence.2", "AA.seq", "Pep.pos.in.AA.seq", "Parent.Entry")]

## This will take some time...

for (cur.pep.phos in 1:nrow(TMT.results)) {
  
  setTxtProgressBar(progress.bar.2, cur.pep.phos)
  
  pep.rows <- TMT.meta[cur.pep.phos, ]
  
  ## We can have several phosphorylation sites per peptide
  
  phos.site.pep <- vector()
  
  for (cur.row in 1:nrow(pep.rows)) {
    
    phos.sites.pep.temp <- unlist(gregexpr("p", pep.rows[cur.row, "Modified.sequence.2"]))
    
    ## "p" changes position of each "p" by 1 in modified sequences, need to correct
    
    if (length(phos.sites.pep.temp) >= 2) {
      
      phos.pos <- seq(2,length(phos.sites.pep.temp), 1)
      
      norm <- (phos.pos - 1)
      norm <- c(0, norm)
      
      phos.sites.pep.temp.2 <- phos.sites.pep.temp - norm
      phos.site.pep <- c(phos.site.pep, phos.sites.pep.temp.2)
      
    } else {
      
      phos.site.pep <- c(phos.site.pep, phos.sites.pep.temp)
      
    }
  }
  
  ## Pep.pos.in.AA.seq -> str_locate gives a list with start and end position, first entry is start position
  phos.sites.abs <- unique(phos.site.pep) + pep.rows$Pep.pos.in.AA.seq[1] - 1
  
  ## Which amino acid is phosphorylated?
  
  if (length(unique(phos.site.pep)) >= 2) {
    
    phos.aa <- vector()
    
    for (e in 1:length(unique(phos.sites.pep.temp))) {
      
      phos.aa.temp <- substring(pep.rows$Modified.sequence.2, phos.sites.pep.temp[e] + 1, phos.sites.pep.temp[e] + 1)
      
      phos.aa <- c(phos.aa, phos.aa.temp)
      
    }
    
  } else {
    
    phos.aa <- substring(pep.rows$Modified.sequence.2, phos.site.pep + 1, phos.site.pep + 1)
    
  }
  
  phos.sites.abs.aa <- paste(phos.aa, phos.sites.abs, "-p", sep = "")
  
  pep.rows[, "Phosphorylation.sites"] <- paste(phos.sites.abs.aa, collapse = ",")
  
  cur.phos <- paste(phos.sites.abs.aa, collapse = ",")
  
  cur.prot <- pep.rows$Parent.Entry
  
  row.psp <- which(phos.site.plus$ACC_ID == cur.prot)
  
  psp.temp <- phos.site.plus[row.psp, "MOD_RSD"]
  
  phos.flag <- logical()
  
  iterations <- length(unlist(gregexpr("-", cur.phos)))
  
  if (iterations >= 2) {
    
    cur.phos <- unlist(strsplit(cur.phos, split = ","))
    
  }
  
  for (iter in 1:iterations) {
    
    if (cur.phos[iter] %in% psp.temp) {
      
      phos.flag.temp <- TRUE
      phos.flag <- c(phos.flag, phos.flag.temp)
      
    } else {
      
      phos.flag.temp.2 <- FALSE
      phos.flag <- c(phos.flag, phos.flag.temp.2)
    }
    
  }
  
  pep.rows[, "Known.Phos.Site.Plus"] <- paste(phos.flag, collapse = ",")
  
  ## Phosphorylation site part of RPXS/T motif?
  
  temp <- unlist(str_match_all(pep.rows[, "Phosphorylation.sites"], "\\d*"))
  temp.2 <- as.integer(unlist(temp[which(temp != "")]))
  
  pot.sites <- length(temp.2)
  
  rpxs <- logical()
  
  for (cur.site in 1:pot.sites) {
    
    rpxs.temp <- substring(pep.rows[, "AA.seq"], temp.2[cur.site] - 3, temp.2[cur.site])
    
    rpxs.temp.2 <- unlist(grep("^RP[A-Z][S,T]", rpxs.temp))
    
    if (length(rpxs.temp.2) >= 1) {
      
      rpxs.temp.3 <- TRUE
      
      rpxs <- c(rpxs, rpxs.temp.3)
      
    } else {
      
      rpxs.temp.3 <- FALSE
      
      rpxs <- c(rpxs, rpxs.temp.3)
      
    }
    
  }
  
  pep.rows[, "Part.of.RPXS.T.motif"] <- paste(rpxs, collapse = ",")
  
  ## Remove AA sequence, Modified.sequence.2, Parent.Entry
  col.rem <- grep("AA.seq|Modified.sequence.2|Parent.Entry", colnames(pep.rows))
  pep.rows <- pep.rows[, -col.rem]
  
  pep.final <- rbind(pep.final, pep.rows)
  
}

TMT.results <- merge(TMT.results, pep.final, by = "Counter")

## GO terms
uniprot.go.clean <- uniprot.go[, c(1, 5:7, 11)]
colnames(uniprot.go.clean)[1] <- "Parent.Entry"
colnames(uniprot.go.clean)[2] <- "Gene.Names.all"

TMT.results.temp <- merge(TMT.results, uniprot.go.clean, by = "Parent.Entry")

## This merge() dropped some rows, re-estate them

dropped <- setdiff(TMT.results$Counter, TMT.results.temp$Counter)

TMT.dropped <- data.frame()

for (cur.drop in dropped) {
  
  gene.temp <- TMT.results[which(TMT.results$Counter == cur.drop), "Gene.names"]
  gene.temp <- gsub(";", "|", gene.temp)
  go.temp <- grep(paste0("\\b", gene.temp, "\\b"), uniprot.go.clean$Gene.Names.all)

  data.temp <- cbind(TMT.results[which(TMT.results$Counter == cur.drop),], uniprot.go.clean[go.temp, 2:5])
  
  TMT.dropped <- rbind(TMT.dropped, data.temp)
  
}

TMT.results <- rbind.data.frame(TMT.results.temp, TMT.dropped, stringsAsFactors = FALSE)

## Fold change 

TMT.results$Fold.change <- ifelse(
  TMT.results$logFC > 0,
  2^(TMT.results$logFC),
  as.numeric(paste0("-", 2^(abs(TMT.results$logFC))))
)

## Write TMT results file ----
## Remove AA sequence, increases file size unnecessarily

message(paste0("\nWriting ", output.res.dir, results.file))

aa.seq.col <- grep("^AA.seq", colnames(TMT.results))

write.csv(
  TMT.results[, -aa.seq.col], 
  paste0(output.res.dir, results.file),
  row.names = FALSE
)

## Clean up TMT results ----
## Filter by adj.p.value <= 0.65 and sort by fold-change

TMT.results.clean <- TMT.results[which(TMT.results$adj.P.Val <= 0.65), ]

TMT.results.clean <- TMT.results.clean[order(TMT.results.clean$Fold.change, decreasing = TRUE), ]

## Write cleaned up TMT results ----
## Remove AA sequence, increases file size unnecessarily

message(paste0("Writing ", output.res.dir, results.clean.file))

aa.seq.col <- grep("^AA.seq", colnames(TMT.results.clean))

write.csv(
  TMT.results.clean[, -aa.seq.col],
  paste0(output.res.dir, results.clean.file),
  row.names = FALSE
)

## Garden Sprinkler Plot

message(paste0("Writing ", output.res.dir, gs.file))

TMT.results$gp.plot <- sqrt(1 - TMT.results$adj.P.Val)
TMT.results.base <- TMT.results[which(TMT.results$gp.plot <= 0.5), ]
TMT.results.poi <- TMT.results[which(TMT.results$gp.plot >= 0.5), ]
TMT.results.poi.pos <- TMT.results.poi[which(TMT.results.poi$Fold.change > 0), ]

## Name of validated proteins

val.pep.names <- c(
  "CDKL5",
  "YLPM1",
  "EP400",
  "MPLKIP",
  "TCEB3",
  "ZNF592"
)

val.pep <- which(TMT.results.poi.pos$Gene.names %in% val.pep.names)

TMT.results.poi.pos <- TMT.results.poi.pos[val.pep, ]

## Adjust gene names to synonyms to harmonize names in manuscript

TMT.results.poi.pos$Gene.names <- sub("YLPM1", "ZAP3", TMT.results.poi.pos$Gene.names)
TMT.results.poi.pos$Gene.names <- sub("MPLKIP", "TTDN1", TMT.results.poi.pos$Gene.names)
TMT.results.poi.pos$Gene.names <- sub("TCEB3", "ELOA", TMT.results.poi.pos$Gene.names)

##Sprinkler plot Figure 2B ----

pal.2 <- rev(wes_palette("Darjeeling1"))
pal.2 <- pal.2[c(1,4,3,2,5)]

ggplot(
  TMT.results,
  aes(x = logFC, y = gp.plot)
  ) +
  geom_pointdensity() +
  scale_color_gradientn(
    colours = pal.2,
    name = "Number of \nneighbouring \npeptides",
    label = comma
  ) + 
  geom_point(
    data = TMT.results.poi.pos,
    colour = "black",
    size = 2
  ) +
  geom_text_repel(
    data = TMT.results.poi.pos,
    aes(label = Gene.names),
    min.segment.length = 0, 
    family = "Arial",
    fontface = "bold",
    size = 2
  ) +
  ylab(expression(sqrt(1-adj.~italic(p)-Value))) +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(log[2]~Fold~Change)) +
  theme(
    axis.title = element_text(size = 12),
    axis.text.x = element_text(family = "Arial", size = 10),
    axis.text.y = element_text(family = "Arial", size = 10),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = c(0.15, 0.77),
  )

ggsave(
  paste0(output.res.dir, gs.file),
  width = 109,
  height = 109,
  unit = "mm",
  dpi = 600
)

## Sprinkler plot without labels

ggplot(
  TMT.results,
  aes(x = logFC, y = gp.plot)
) +
  geom_pointdensity() +
  scale_color_gradientn(
    colours = pal.2,
    name = "Number of \nneighbouring \npeptides",
    label = comma
  ) + 
  geom_point(
    data = TMT.results.poi.pos,
    colour = "black",
    size = 2
  ) +
  ylab("") +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(log[2]~Fold~Change)) +
  theme(
    axis.title = element_text(size = 12),
    axis.text.x = element_text(family = "Arial", size = 10),
    axis.text.y = element_text(family = "Arial", size = 10),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = c(0.15, 0.77),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

ggsave(
  paste0(output.res.dir, gs.file.2),
  width = 109,
  height = 109,
  unit = "mm",
  dpi = 1200,
  bg = "transparent"
)


### Descriptive statistics dataframe ----

message(paste0("Writing ", output.res.dir, stats.file))

stats.df <- data.frame(
  Group = character(),
  Number = numeric()
)

total.id <- nrow(data.mq[-c(rem.rev, rem.cont, rem.pep ), ])
total.id.pif <- nrow(data.clean)
ratio.pif <- round((total.id.pif / total.id), 4) * 100
ratio.incomplete.kd.data <- round(1 - (nrow(data.clean) / nrow(data.preclean)), 4) * 100
pep.data <- nrow(data.preclean)
pep.quant <- nrow(data.av)
non.phos <- length(rem.non.phos) 

ratio.cont <- round(
  length(unique(data.av[rem.non.phos, "Modified.sequence"])) /  
  length(unique(data.av$Modified.sequence)),
  4
) * 100

data.quant <- nrow(data.limma)
uni.seq.quant <- length(unique(TMT.results$Modified.sequence.2))

uni.phos.seq <- TMT.results[!duplicated(TMT.results$Modified.sequence.2), ]

no.phos.sites <- sum(as.numeric(uni.phos.seq$Phospho.STY)) 

stats.df <- data.frame(
  Group = c("Total number of peptides identified (5%FDR)",
            "Total number of peptides with PIF >= 0.60",
            "Percentage of peptides with PIF >= 0.60",
            "Percentage of peptides with < 2 Observations in KD|WT Group",
            "Number of peptides before averaging",
            "Number of averaged peptides quantified",
            "Number of non-phospho peptides",
            "Percentage of non-phospho peptides",
            "Number of phospho peptides statistically tested",
            "Number of unique phospho peptide sequences quantified",
            "Number of unique phosphorylation sites quantified"
  ),
  Number = c(round(total.id, 1),
             total.id.pif,
             ratio.pif,
             ratio.incomplete.kd.data,
             pep.data,
             pep.quant,
             non.phos,
             ratio.cont,
             data.quant,
             uni.seq.quant,
             no.phos.sites
  )
)

write.csv(
  stats.df,
  paste0(output.res.dir, stats.file),
  row.names = FALSE
)

## Table Figure 2C ----

cols.fig <- c(
  "Gene.names",
  "protein.name.long",
  "Parent.Entry",
  "Fold.change",
  "Phosphorylation.sites",
  "Phospho.STY.probabilities",
  "Modified.sequence.2"
)

table.fig <- TMT.results.clean[which(TMT.results.clean$Fold.change > 0), cols.fig]
table.fig$Fold.change <- round(table.fig$Fold.change, 2)

## Adjust gene and protein names to harmonize naming throughout manuscript

table.fig$Gene.names <- sub("YLPM1", "ZAP3", table.fig$Gene.names)
table.fig$Gene.names <- sub("MPLKIP", "TTDN1", table.fig$Gene.names)
table.fig$Gene.names <- sub("TCEB3", "ELOA", table.fig$Gene.names)

table.fig$protein.name.long <- sub(
  "YLP motif-containing protein 1",
  "Nuclear protein ZAP 3",
  table.fig$protein.name.long
)

table.fig$protein.name.long <- sub(
  "M-phase-specific PLK1-interacting protein",
  "Tricothiodystrophy non-photosensitive 1 protein",
  table.fig$protein.name.long
)

table.fig$protein.name.long <- sub(
  "M-phase-specific PLK1-interacting protein",
  "Tricothiodystrophy non-photosensitive 1 protein",
  table.fig$protein.name.long
)

## Extract phosphorylation motif

table.fig$motif <- substring(
  table.fig$Modified.sequence.2,
  regexpr("p", table.fig$Modified.sequence.2) - 3,
  regexpr("p", table.fig$Modified.sequence.2) + 2
) 

## Clean-up and write table

colnames(table.fig) <- c(
  "Protein",
  "",
  "Uniprot.Acc.No.",
  "Fold.Change",
  "Phosphorylation.Sites",
  "Phospho.STY.probabilities",
  "Peptide.Sequence",
  "Motif"
)

write.csv(
  table.fig,
  paste0(output.res.dir, "Table_Figure_2C.csv"),
  row.names = FALSE
)

message("Analysis done \n")
