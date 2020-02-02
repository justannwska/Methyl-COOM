## 02_Linear modeling on all B-cell specific CpGs.R
## 
## 
##
## 
########################################################################################################################





#########################################################
# libraries
#
library(ape)
library(openxlsx)
library(phylobase)
library(gdata)






###########################################################
# directories
#
analysis.dir <-  "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/analysis/" # output directory for the analysis
data.dir <-  "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/data/" # directory containing linear B cell-specific CpGs

annotation.dir <- "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/annotation/"
phylogeny.dir <- "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/analysis/01_phylogeny/"





###########################################################
# load annotation file for CLL
annotation.file <- read.xlsx(paste0(annotation.dir, "/all_450k_sample_NG_2019-10-31.xlsx"),1)




############################################################
# methylation data
meth.data.all <- readRDS(file.path(data.dir, "avg.meth.bcell.CLL.samples.on.Bcells_specific.probes.RDS"))
meth.data.full <- meth.data.all




## prepare a data frame for modeling
meth.data.without.cg <- meth.data.all
rownames(meth.data.without.cg) <- meth.data.without.cg[ , "cg"]
meth.data.without.cg <- meth.data.without.cg[ ,-which(colnames(meth.data.without.cg)=="cg")]




########################################################
# load the distances for B cells and CLLs determined using package "ape"
Distances.Bcells <- readRDS(paste0(phylogeny.dir, "/Manhattan.phylogeny.normal.B.cells.on.linear.Bcell.specifics.CpGs2020-02-01.RDS"))
Distances.CLL.all <- readRDS(paste0(phylogeny.dir, "/COO.assignment.Chris.NGcohort.on.linear.Bcell.specifics.CpGs2020-02-01.RDS"))

Distances.CLL.all$class <- NA
for(patients in 1:nrow(Distances.CLL.all)){
  Distances.CLL.all$class[patients] <- annotation.file$cluster[Distances.CLL.all$sample[patients] == annotation.file$ALT.NAME]
  
}




########################################################
# calculate regression model for each cpg in meth table
meth.data.without.cg$intmbcs.mgzs <- rowMeans(meth.data.without.cg[, c("intMBCs", "MGZs")])
Meth.Bcells <- meth.data.without.cg[,c("NBCs", "GCFs", "loMBCs", "intmbcs.mgzs", "hiMBCs")]


#Bcells
linear.model <- data.frame(par.a=numeric(nrow(Meth.Bcells)),
                           par.b=numeric(nrow(Meth.Bcells))
)
for(rows in 1:nrow(Meth.Bcells)){
  data.for.regression <- cbind.data.frame(t(Meth.Bcells[rows,]), Distances.Bcells$Location)
  colnames(data.for.regression) <- c("Meth", "Distances")
  data.for.regression$Meth <- data.for.regression$Meth *100
  linear.regression <- lm(Meth ~ Distances, data.for.regression)
  linear.model$par.a[rows] <- linear.regression$coefficients[2]
  linear.model$par.b[rows] <- linear.regression$coefficients[1]
  
}
rownames(linear.model) <- rownames(Meth.Bcells)
saveRDS(linear.model, paste0(data.dir, "/linear.model.on.B.cell.specific.CpGs.Bcells-", Sys.Date(), ".RDS"))










# linear model for CLLs 
Meth.CLLs.all <- meth.data.without.cg[ , Distances.CLL.all$sample]
Matrix.for.relative.meth.CLL.group1 <- list()


for(rows in 1:nrow(Meth.CLLs.all)){
  
  cg.id <- rownames(Meth.CLLs.all)[rows]
  Matrix.for.relative.meth.CLL.group1[[rows]] <- data.frame(Meth=numeric(ncol(Meth.CLLs.all)),
                                                            Distances=numeric(nrow(Distances.CLL.all)),
                                                            Meth.cell.of.origin = numeric(nrow(Distances.CLL.all)),
                                                            Relative.meth =numeric(nrow(Distances.CLL.all)))
  
  
  Matrix.for.relative.meth.CLL.group1[[rows]]$Meth <- c(t(Meth.CLLs.all[rows,]))
  names <- rownames(t(Meth.CLLs.all[rows,]))
  rownames(Matrix.for.relative.meth.CLL.group1[[rows]]) <- names
  
  Matrix.for.relative.meth.CLL.group1[[rows]]$Meth <- Matrix.for.relative.meth.CLL.group1[[rows]]$Meth *100
  Matrix.for.relative.meth.CLL.group1[[rows]]$Distances <- Distances.CLL.all$location
  bcell.regression.model <- linear.model[rownames(linear.model) == cg.id, ]
  Matrix.for.relative.meth.CLL.group1[[rows]]$Meth.cell.of.origin <- Matrix.for.relative.meth.CLL.group1[[rows]]$Distances*bcell.regression.model$par.a + bcell.regression.model$par.b
  Matrix.for.relative.meth.CLL.group1[[rows]]$Relative.meth <- Matrix.for.relative.meth.CLL.group1[[rows]]$Meth - Matrix.for.relative.meth.CLL.group1[[rows]]$Meth.cell.of.origin
  colnames(Matrix.for.relative.meth.CLL.group1[[rows]]) <- c("Meth","Distances", "Meth.cell.of.origin","Relative.meth")

}
names(Matrix.for.relative.meth.CLL.group1) <- rownames(Meth.CLLs.all)




## Determine mean relative methylation for all CLLs
Mean.relative.meth.CLL.all <- lapply(Matrix.for.relative.meth.CLL.group1, function(get.mean.values){
  mean(get.mean.values$Relative.meth)
})

Mean.meth.CLL.all <- do.call(rbind.data.frame, Mean.relative.meth.CLL.all)
colnames(Mean.meth.CLL.all) <- "Relative Methylation"
rownames(Mean.meth.CLL.all) <- names(Matrix.for.relative.meth.CLL.group1)
  
  
saveRDS(Matrix.for.relative.meth.CLL.group1, paste0(data.dir, "/linear.model.on.B.cell.specific.CpGs.CLLs-",Sys.Date(), ".RDS"))
saveRDS(Mean.relative.meth.CLL.all, paste0(data.dir, "/mean.meth.linear.model.on.B.cell.specific.CpGs.CLLs-",Sys.Date(), ".RDS"))







## Prepare relative methylation events for CLLs
Relative.meth.per.cg.all.patients  <- lapply(Matrix.for.relative.meth.CLL.group1, function(get.meth){
  data.frame(t(get.meth$Relative.meth))
})
rel.meth <- do.call(rbind.data.frame, Relative.meth.per.cg.all.patients)
colnames(rel.meth) <- rownames(Matrix.for.relative.meth.CLL.group1[[1]])

na.values <- apply(rel.meth, 1, function(x){any(is.na(x))})
na.values.df <- rel.meth[na.values, ]
rel.meth.non.zeros <- rel.meth[!na.values, ]



## save relative methylation calls into .RDS file
saveRDS(rel.meth, paste0(data.dir, "/Relative.meth.linear.model.on.B.cell.specific.CpGs.CLLs-",Sys.Date(), ".RDS"))









