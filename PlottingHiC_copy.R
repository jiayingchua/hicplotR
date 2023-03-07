## Load packages
#library("org.Hs.eg.db")
#library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#library("plotgardenerData")
#library("AnnotationHub")
#library(plotgardenerData)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")

BiocManager::install("plotgardener")
library("plotgardener")
library(RColorBrewer)

## Load files
hiCFile <-
  "C://Users//JiaYing//Group Project//Ridaeus_Ras1_scaffolds_yahs_JBAT.hic"
assemblyfile <-
  "C://Users//JiaYing//Group Project//Ridaeus_Ras1_scaffolds_yahs_JBAT.curated4.assembly"


######### Checks contents of hic file, converts hic into df ########

chromname <- strawr::readHicChroms(hiCFile)
bpreso <- strawr::readHicBpResolutions(hiCFile)
hiCdf <-
  strawr::straw("NONE", hiCFile, "assembly", "assembly", "BP", 250000, "observed")
# plots entire assembly
x11()
hiCPlot <-
  plotgardener::plotHicSquare(
    data = hiCdf,
    chrom = "chr1",
    norm = "NONE",
    palette = colorRampPalette(brewer.pal(n = 9, "Reds")),
    colorTrans = "log10"
  )

######### Convert assembly file into appropriate lists/vectors #########

processFile = function(filepath) {
  contigs = c()
  contignum = c()
  pixels = c()
  chromosomes = c()
  con = file(filepath, "r")
  while (TRUE) {
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    if (startsWith(line, ">")) {
      #print(paste(line, "TRUE"))
      contigs <- append(contigs, unlist(strsplit(line, " "))[1])
      contignum <-
        append(contignum, unlist(strsplit(line, " "))[2])
      pixels <- append(pixels, unlist(strsplit(line, " "))[3])
    }
    else {
      chromosomes <- append(chromosomes, line)
    }
  }
  assembly <- data.frame(contigs, contignum, pixels)
  return(list(assembly, chromosomes))
  close(con)
}

assembly_chromosomes <- processFile(assemblyfile)

sum_bp <- c(1)
for (i in 1:length(assembly_chromosomes[[2]])) {
  chr <- unlist(strsplit(assembly_chromosomes[[2]][i], " "))
  chr_bp <- c()
  for (j in 1:length(chr)) {
    ct <- abs(as.numeric(chr[j]))
    bp <- as.numeric(assembly_chromosomes[[1]][ct, 3])
    chr_bp <- append(chr_bp, bp)
  }
  sum_bp <- append(sum_bp, sum(chr_bp))
}
print(sum_bp)

## Example: User wants to display chromosome 3
c = 3
chr <- unlist(strsplit(assembly_chromosomes[[2]][c], " "))
chr_bp <- c()
for (j in 1:length(chr)) {
  ct <- abs(as.numeric(chr[j]))
  bp <- as.numeric(assembly_chromosomes[[1]][ct, 3])
  chr_bp <- append(chr_bp, bp)
}
sum_chr <- sum(chr_bp)
print(sum_chr)

########## Plot contact map with scaffold boxes - selected chromosome only ############
hiCdf <-
  strawr::straw("NONE", hiCFile, "assembly", "assembly", "BP", 50000, "observed")

x11()
pageCreate(
  width = 5,
  height = 5,
  default.units = "inches",
  showGuides = FALSE,
  xgrid = 0,
  ygrid = 0
)
hicPlot <- plotHicSquare(
  data = hiCdf,
  chrom = paste("chr", c),
  chromstart = sum(sum_bp[1:c]) - 1,
  chromend = sum(sum_bp[1:c + 1]),
  x = 0,
  y = 0,
  width = 5,
  height = 5,
  default.units = "inches",
  palette = colorRampPalette(brewer.pal(n = 9, "Reds")),
  colorTrans = "log10"
)
annoGenomeLabel(plot = hicPlot,
                scale = "bp",
                x = 0,
                y = 5)

#chr_bp
#sum_chr
newx = 0
newy = 5
for (i in 1:length(chr_bp)) {
  x = chr_bp[i] * (5 / sum_chr)
  y = chr_bp[i] * (5 / sum_chr)
  width = chr_bp[i] * (5 / sum_chr)
  height = -chr_bp[i] * (5 / sum_chr)
  plotRect(
    x = newx,
    y = newy,
    width = width,
    height = height,
    just = c("left", "top"),
    default.units = "inches",
    lwd = 2,
    fill = NA
  )
  newx = newx + x
  newy = newy - y
}

