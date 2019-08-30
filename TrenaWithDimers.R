#adapted from trenaSGM test_NoDnaModelBuilder.R

library(RUnit)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
  load(system.file(package="trenaSGM", "extdata", "mayo.tcx.new.RData"))
#------------------------------------------------------------------------------------------------------------------------

#define functions for construction of a new TF DIMER row from two existing TF (monomer) rows

#takes the minimum of the two rows #passes all tests
minRow <- function(row1, row2) 
{ 
  result <- row1
  for(i in 1:length(row1)) {
    if(row2[i]<row1[i]) {
      result[i] <- row2[i]
    }
  }
  return(result)
}

#takes the square root of the dot product of the two rows #does not successfully build a model
sqrtDotRow <- function(row1, row2) 
{
  result <- row1
  for(i in 1:length(row1)) {
    result[i] <- sqrt(sqrt(abs(row1[i])) * sqrt(abs(row2[i])))
  }
  return(result)
}

#takes the dot product of the two rows #does not succesfully build a model
dotRow <- function(row1, row2) 
{
  result <- sqrt(abs(new.mtx[tf1,]) * abs(new.mtx[tf2,])) #dot product
  return(result)
}

#a test function to mess around with
testRow <- function (row1, row2) { 
  result <- row1
  for(i in 1:length(row1)) {
    if(i > length(row1)/2) {
      result[i] <- row2[i]
    }
  }
  return(result)
}

#------------------------------------------------------------------------------------------------------------------------
#append dimers to mtx

#read in the expression matrix and list of all known TFs so that we can alter them to include DIMERS
new.mtx <- mtx
new.allTFs <- allKnownTFs()

#read in table of known TF DIMERS
# dimers <- read.table(file = '/Users/user/Downloads/TF_DIMERS_PREPPED.tsv', sep = '\t', header = TRUE)

#append a fake tf dimer to fit in the matrix
tf1 <- "IRF5"
tf2 <- "IKZF1"
newRow <- minRow(new.mtx[tf1,], new.mtx[tf2,])
new.mtx <- rbind(newRow, new.mtx)
dimername <- paste(tf1, "_", tf2, sep="")
rownames(new.mtx)<-c(dimername ,rownames(new.mtx)[-1])
new.allTFs <- c(dimername, new.allTFs)
test.mtx <- new.mtx[c("IRF5","IKZF1", "IRF5_IKZF1"),] #to check and make sure the algorithm worked

#loop through all known TF DIMERS and append new rows to the mtx
# for (i in 1:nrow(dimers)) {
#   tf1 <- as.character(dimers[i, 1])
#   tf2 <- as.character(dimers[i,2])
#   newRow <- minRow(new.mtx[tf1,], new.mtx[tf2,])
#   new.mtx <- rbind(newRow, new.mtx)
#   dimername <- paste(tf1, "_", tf2, sep="")
#   rownames(new.mtx)<-c(dimername ,rownames(new.mtx)[-1])
#   new.allTFs <- c(dimername, new.allTFs)
# }

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_constructor()
  test_build.trem2.noDNA.13.known.TFs()
  test_build.trem2.noDNA.all.known.TFs()
  test_build.trem2.noDNA.bogus.targetGene()
} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
  printf("--- test_constructor")
  
  genome <- "hg38"
  targetGene <- "TREM2"
  chromosome <- "chr6"
  tss <- 41163186
  # strand-aware start and end: trem2 is on the minus strand
  tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)
  
  build.spec <- list(title="trem2.rmm.2000up.200down",
                     type="noDNA.tfsSupplied",
                     candidateTFs=c("HLF", "STAT4", "SATB2", "SATB1", "TSHZ3", "TSHZ2", "FOXP2", "IRF5_IKZF1"),
                     tfPool=new.allTFs,#altered version of allKnownTFs()
                     matrix=new.mtx, #altered version of mtx, as loaded from mayo.tcx.new.RData
                     tfPrefilterCorrelation=0.4,
                     annotationDbFile=dbfile(org.Hs.eg.db),
                     orderModelByColumn="rfScore",
                     solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                     quiet=TRUE)
  
  builder <- NoDnaModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
  
  checkTrue("NoDnaModelBuilder" %in% is(builder))
  
} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.noDNA.13.known.TFs <- function()
{
  printf("--- test_build.trem2.noDNA.13.known.TFs")
  
  genome <- "hg38"
  targetGene <- "TREM2"
  
  candidate.tfs <- c("IRF5", "IKZF1", "LYL1", "SPI1", "CEBPA", "TFEC",
                     "BHLHE41", "IRF8", "TAL1","ELK3", "POU2F2", "MAFB",
                     "ZBTB18", "bogus", "IRF5_IKZF1") #altered to include the made up dimer IRF5_IKZF1 as 
  
  build.spec <- list(title="trem2.noDNA.allTFs",
                     type="noDNA.tfsSupplied",
                     matrix=new.mtx, #mtx with dimer rows appended
                     candidateTFs=candidate.tfs, 
                     tfPool=new.allTFs, #altered version of allKnownTFs()
                     tfPrefilterCorrelation=0.2,
                     annotationDbFile=dbfile(org.Hs.eg.db),
                     orderModelByColumn="rfScore",
                     solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                     quiet=TRUE)
  
  builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
  x <- build(builder)
  checkEquals(x$regulatoryRegions, data.frame())
  tbl.model <- x$model
  # with a relaxed tfPrefilterCorrelation, and the hand-picked TFs listed above
  # all but "bogus" make the cut
  checkEquals(setdiff(candidate.tfs, tbl.model$gene), "bogus")
  # noDNA implies no bindingSites
  checkTrue(all(is.na(tbl.model$bindingSites)))
  
  # all pearsonCoeff above prefilter threshold?
  checkTrue(all(abs(tbl.model$pearsonCoeff) > 0.2))
  
  
  #--------------------------------------------------------------------------------
  # run again with a stricter tfPrefilterCorrelation
  #--------------------------------------------------------------------------------
  build.spec$tfPrefilterCorrelation <- 0.7
  builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
  x <- build(builder)
  tbl.model <- x$model
  checkTrue(all(abs(tbl.model$pearsonCoeff) > 0.7))
  checkTrue(all(tbl.model$gene %in% candidate.tfs))
  
} # test_build.trem2.noDNA.13.known.TFS
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.noDNA.all.known.TFs <- function()
{
  printf("--- test_build.trem2.noDNA.all.known.TFs")
  
  genome <- "hg38"
  targetGene <- "TREM2"
  
  candidate.tfs <- new.allTFs #altered version of allKnownTFs()
  
  build.spec <- list(title="trem2.noDNA.allTFs",
                     type="noDNA.tfsSupplied",
                     matrix=new.mtx, #mtx with dimer rows appended
                     candidateTFs=candidate.tfs,
                     tfPool=new.allTFs, #altered version of allKnownTFs()
                     tfPrefilterCorrelation=0.7,
                     annotationDbFile=dbfile(org.Hs.eg.db),
                     orderModelByColumn="pearsonCoeff",
                     solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                     quiet=TRUE)
  
  builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
  x <- build(builder)
  tbl.model <- x$model
  checkEquals(dim(x$regulatoryRegions), c(0,0))
  checkTrue(all(tbl.model$peasonCoeff > 0.7))
  # the order
  checkEquals(tbl.model$gene, c("IRF5_IKZF1", "PLEK", "IRF5", "IKZF1", "LYL1", "SPI1", "TFEC"))
  
} # test_build.trem2.noDNA.all.known.TFs
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.noDNA.bogus.targetGene <- function()
{
  printf("--- test_build.trem2.noDNA.bogus.targetGene")
  
  genome <- "hg38"
  targetGene <- "bogusGene"
  
  candidate.tfs <- new.allTFs #altered version of allKnownTFs()
  
  build.spec <- list(title="trem2.noDNA.allTFs",
                     type="noDNA.tfsSupplied",
                     matrix=new.mtx, #mtx with dimer rows appended
                     candidateTFs=candidate.tfs,#altered version of allKnownTFs()
                     tfPool=new.allTFs, #altered version of allKnownTFs()
                     annotationDbFile=dbfile(org.Hs.eg.db),
                     tfPrefilterCorrelation=0.7,
                     orderModelByColumn="pearsonCoeff",
                     solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                     quiet=TRUE)
  
  checkException(builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE), silent=TRUE)
  
} # test_build.trem2.noDNA.bogus.targetGene
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()
