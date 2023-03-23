## PLOTS ALL 7 CHROMOSOMES AND SAVES EACH CHR TO ONE BMP/PDF/SVG FILE
## INPUTS - .hic file and .assembly file, resolution, noc

# rm(list = ls())
# graphics.off()

args = commandArgs(trailingOnly = TRUE)
# argument 1 = hic file
# argument 2 = assembly file
# argument 3 = output folder path
# argument 4 = resolution
# argument 5 = number of chromosomes
hiCFile <- args[1]
assemblyfile <- args[2]
output_folder <- args[3]
reso <- as.numeric(args[4])
noc <- as.numeric(args[5])


## Load packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# BiocManager::install("plotgardener")
# install.packages("strawr")
# library("plotgardener")
# library(RColorBrewer)

# List of packages for session
.packages = c("plotgardener", "RColorBrewer")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

## Load files
# hiCFile <-
#   "C://Users//JiaYing//Group Project//Ridaeus_Ras1_scaffolds_yahs_JBAT.hic"
# assemblyfile <-
#   "C://Users//JiaYing//Group Project//Ridaeus_Ras1_scaffolds_yahs_JBAT.curated4.assembly"
# assemblyfile <- 
#   "C://Users//JiaYing//Group Project//Ridaeus_Ras1_v1.0//Ridaeus_Ras1_scaffolds_yahs_JBAT.assembly"

# Specify output folder
# output_folder <- "C://Users//JiaYing//Group Project//HiC_images"
dir.create(output_folder)


# number of chromosomes (noc)
#noc <- 7

######### Checks contents of hic file ########

chromname <- strawr::readHicChroms(hiCFile)
bpreso <- strawr::readHicBpResolutions(hiCFile)
#reso <- 50000 #resolution (user input)
width <- 5 #inches
height <- 5 #inches

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

# loop for number of chromosomes

sum_bp <- c(1)
for (i in 1:noc) {
  chr <- unlist(strsplit(assembly_chromosomes[[2]][i], " "))
  chr_bp <- c()
  for (j in 1:length(chr)) {
    ct <- abs(as.numeric(chr[j]))
    bp <- as.numeric(assembly_chromosomes[[1]][ct, 3])
    chr_bp <- append(chr_bp, bp)
  }
  sum_bp <- append(sum_bp, sum(chr_bp))
}

## Plot for all chromosomes
for (c in 1:noc) {
  chr <- unlist(strsplit(assembly_chromosomes[[2]][c], " "))
  chr_bp <- c()
  for (j in 1:length(chr)) {
    ct <- abs(as.numeric(chr[j]))
    bp <- as.numeric(assembly_chromosomes[[1]][ct, 3])
    chr_bp <- append(chr_bp, bp)
  }
  sum_chr <- sum(chr_bp)
  
  ## Plot contact map with scaffold boxes - selected chromosome only
  hiCdf <-
    strawr::straw("NONE",
                  hiCFile,
                  "assembly",
                  "assembly",
                  "BP",
                  reso,
                  "observed")
  
  png(paste(output_folder,"/hicplot_chr", c, ".png", sep=""), # File name
    width = width,
    height = height,
    units = "in", # required for .bmp and .png
    res = 100, # required for .bmp and .png
    bg = "white" # Background color
  )         
  
  pageCreate(
    width = width,
    height = height,
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
    width = width,
    height = height,
    default.units = "inches",
    palette = colorRampPalette(brewer.pal(n = 9, "Reds")),
    colorTrans = "log10"
  )
  # annoGenomeLabel(
  #   plot = hicPlot,
  #   scale = "bp",
  #   x = 0,
  #   y = height
  # )
  
  #chr_bp
  #sum_chr
  newx = 0
  newy = height
  for (i in 1:length(chr_bp)) {
    x = chr_bp[i] * (width / sum_chr)
    y = chr_bp[i] * (height / sum_chr)
    x_width = chr_bp[i] * (width / sum_chr)
    y_height = -chr_bp[i] * (height / sum_chr)
    plotRect(
      x = newx,
      y = newy,
      width = x_width,
      height = y_height,
      just = c("left", "top"),
      default.units = "inches",
      lwd = 1,
      fill = NA
    )
    newx = newx + x
    newy = newy - y
  }
  
  dev.off()
}
