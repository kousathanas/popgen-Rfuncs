####################################
#read in plink genotype file
####################################
read.plink.output.genotypeM<-function(file){
  require(snpStats)
  require(BEDMatrix)
  
  pathM <- paste(file, c(".bed", ".bim", ".fam"), sep = "")
  SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
  x<-as.matrix(BEDMatrix(paste(file,".bed",sep='')))
  return(x)
}
