## 05_Determine class C & D events.R
## 
## 
## 
## These events are calculated on non B cell specific CpGs
########################################################################################################################






#########################################################
# libraries
#
library(openxlsx)
library(gdata)
#library(ggplot2)
library(ggpubr)
library(reshape2)




###########################################################
# directories
#
analysis.dir <- "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/analysis/"
data.dir <-  "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/data/"
annotation.dir <- "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/annotations/"
phylogeny.dir <-  "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/analysis/01_phylogeny/"
class.dir <- "/Users/justynaanna/Documents/CLL/Github/MethylCOOM/analysis/"







######################################################################
# load relative methylation changes based on calculated cell-of-origin
#
relative.methylation.changes <- readRDS(paste0(linear.dir, "/Relative.meth.linear.model.on.non.B.cell.specific.CpGs.CLLs-2020-01-03.RDS"))






#######    Determine aberrant disease-specific methylation events (with a cutoff of 20%)  ####################################
# frequencies methylation change more than 20%  methylation events (future class D)
frequencies.true <- data.frame("nr.patients" = as.numeric())
for(element in 1:nrow(relative.methylation.changes)){
  print(element)
  table.frequencies <- table(relative.methylation.changes[element, ] > 20 ) 
  
  if(any(names(table.frequencies) )== "TRUE"){
    frequencies.true[element, "nr.patients"] <- table.frequencies[["TRUE"]]
  } else{ 
    print(element)
    frequencies.true[element, "nr.patients" ] <- 0
  }
}

rownames(frequencies.true) <-  rownames(relative.methylation.changes)
frequencies.true$cg <- rownames(frequencies.true)



# optional: plot histogram for methylation changes
# plot(density(frequencies.true$nr.patients))
# hist(frequencies.true$nr.patients, xlab= "Number of patients", ylab="Number of CpG sites", main="", 
#      cex=1, cex.lab=1, ylim=c(0,10000), las=2)
# abline(v=31, col="blue", lwd=2)

#sample.fraction <- quantile(frequencies.true$nr.patients)




# frequencies methylation change less than 20% (future class C)
frequencies.true.less <- data.frame("nr.patients" = as.numeric())
for(element in 1:nrow(relative.methylation.changes)){
  print(element)
  table.frequencies <- table(relative.methylation.changes[element, ] < (-20) ) 
  
  if(any(names(table.frequencies) )== "TRUE"){
    frequencies.true.less[element, "nr.patients"] <- table.frequencies[["TRUE"]]
  } else{ 
    print(element)
    frequencies.true.less[element, "nr.patients" ] <- 0
  }
}

rownames(frequencies.true.less) <-  rownames(relative.methylation.changes)
frequencies.true.less$cg <- rownames(frequencies.true.less)




#optional: plotting histogram for methylation changes
# plot(density(frequencies.true.less$nr.patients))
# hist(frequencies.true.less$nr.patients, xlab= "Number of patients", ylab="Number of CpG sites", main="", 
#      cex=1, cex.lab=1, ylim=c(0,10000), las=2)
# abline(v=31, col="blue", lwd=2)
# 
# sample.fraction <- quantile(frequencies.true.less$nr.patients)
# x <- names(sample.fraction)
# plot(x ~ sample.fraction)
# 





###################################################################
# save the files with frequencies
write.csv(frequencies.true, paste0(class.dir, "//frequencies.classD.more.than.0.2 cutoff.aberrant.methylation.events.per.patient", Sys.Date(), ".csv"))
write.csv(frequencies.true.less, paste0(class.dir, "/frequencies.classC.less.than.0.2 cutoff.aberrant.methylation.events.per.patient", Sys.Date(), ".csv"))






####################################################################################
# select the CpG sites with 20%  cutoff for variable number of patients (90% - 10% of patients)
# Later on one can decide which cutoff for the number of patients to use
number.of.patients <- ncol(relative.methylation.changes)

## determine an absolute number of patients belonging to each threshold
patient.fraction <- c( 0.9, 0.8, 0.75, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1)
patient.subset <- round(number.of.patients * patient.fraction)


## looping over all frequencies
classC.list <- list()
classD.list <- list()

relative.methylation.classC <- list()
relative.methylation.classD <- list()



for(element in 1:length(patient.subset)) {
  
  element.no <- patient.subset[[element]]
  
  ## determine number of CpGs in classA & B
  classC.list[[element]] <- data.frame(frequencies.true.less[frequencies.true.less$nr.patients >= as.numeric(element.no), ])
  classC.list[[element]]$cg <- as.character(classC.list[[element]]$cg)
  
  classD.list[[element]] <- data.frame(frequencies.true[frequencies.true$nr.patients >= as.numeric(element.no), ])
  classD.list[[element]]$cg <- as.character(classD.list[[element]]$cg)
  
  
  
  ## extract relative methylation values for the classC & E events
  
  # classC
  relmeth.objectD <- relative.methylation.changes[rownames(relative.methylation.changes) %in% classC.list[[element]]$cg, ]
  relmeth.objectD$cpg <- rownames(relmeth.objectD)
  
  if( all(relmeth.objectD$cpg %in% classC.list[[element]]$cg) ) {
    print("")
  } else {
    print("error: relative methylation was not extracted for all CpGs from classD")
  }
  
  # barplot per patient all CpGs
  data.for.barplot <- melt(relmeth.objectD[1:4,], id="cpg")
  
  p = ggplot(data = data.for.barplot, aes(x = cpg, y = value))
  p = p + geom_bar(stat='identity', width = 0.1) + geom_hline(yintercept= c(-20, 20), color="red")
  p = p + facet_grid(~variable)   
  p
  
  ggsave(paste0(class.dir, "/classC.perCpG.barplot", patient.fraction[element], Sys.Date(), ".pdf"),  plot= p, width= 30)
  
  
  # classD
  relmeth.objectE <- relative.methylation.changes[rownames(relative.methylation.changes) %in% classD.list[[element]]$cg, ]
  relmeth.objectE$cpg <- rownames(relmeth.objectE)
  
  if( all(relmeth.objectD$cpg %in% classD.list[[element]]$cg) ) {
    print("")
  } else {
    print("error: relative methylation was not extracted for all CpGs from classD")
  }
  
  # barplot per patient all CpGs
  data.for.barplot <- melt(relmeth.objectE[1:4,], id="cpg")
  
  p = ggplot(data = data.for.barplot, aes(x = cpg, y = value))
  p = p + geom_bar(stat='identity', width = 0.1) + geom_hline(yintercept= c(-20, 20), color="red")
  p = p + facet_grid(~variable)   
  p
  
  ggsave(paste0(class.dir, "/classD.perCpG.barplot", patient.fraction[element], Sys.Date(), ".pdf"),  plot= p, width= 30)
  
  # assign relative methylation to the list
  relative.methylation.classC[[element]] <- relmeth.objectD
  relative.methylation.classD[[element]] <- relmeth.objectE
  
  
  #check for duplicated CpGs (if duplicated an error will occur)
  #classA
  if(any(duplicated( relative.methylation.classC[[element]]$cpg)) )  {
    print("Error: Duplicated CpGs in extracted relative methylaion classC")
  } else {
    print("")
  }
  
  #classB
  if(any(duplicated(relative.methylation.classD[[element]]$cpg)) ) {
    print("Error: Duplicated CpGs in extracted relative methylaion classD")
  } else {
    print("")
  }
  
  
}

names(relative.methylation.classC) <- c( 0.9, 0.8, 0.75, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1)
names(relative.methylation.classD) <- c( 0.9, 0.8, 0.75, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1)



## saving a file with relative methylation for CpGs defined in classC and D
saveRDS(relative.methylation.classC, paste0(class.dir, "/classC_relative_methylation-", Sys.Date(), ".RDS"))
saveRDS(relative.methylation.classD, paste0(class.dir, "/classD_relative_methylation-", Sys.Date(), ".RDS"))





## extracting number of CLL-specific CpGs, depending on the threshold used 
CpG.count.classC <- sapply(relative.methylation.classC, function(get.cpg.count) { nrow(get.cpg.count) } )
CpG.count.classD <- sapply(relative.methylation.classD, function(get.cpg.count) { nrow(get.cpg.count) } )




## Creating a summary table with a number of events per class per number of patients
summary.table <- data.frame(patient.fraction*100, patient.subset, CpG.count.classC, CpG.count.classD)
colnames(summary.table) <- c("Threshold", "No_Patients", "No_CpGs_classC", "No_CpGs_classD")
summary.table[, c("No_CpGs_classC", "No_CpGs_classD")] <- apply(summary.table[, c("No_CpGs_classC", "No_CpGs_classD")], 2, function(x) as.numeric(as.character(x)));
summary.table$Threshold <- as.character(summary.table$Threshold)
write.xlsx(summary.table, paste0(class.dir, "/classC&E_statistics", Sys.Date(), ".xlsx"))




## plot the summary results
summary.data <- data.frame(rep(patient.fraction*100,2), c(CpG.count.classC, CpG.count.classD))
summary.data$class <- c( rep("classC", length(patient.subset)), rep("classD", length(patient.subset)) )
colnames(summary.data) <- c("Threshold", "No_CpGs", "Class")
summary.data$Threshold <- as.numeric(summary.data$Threshold)



# barplot
pdf(paste0(class.dir,"classC.statistics.barplot", Sys.Date(), ".pdf"), height=4, width=5)
ggbarplot(summary.table, x = "Threshold", 
          y = "No_CpGs_classC",
          fill = "dodgerblue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          #palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          #sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "No. CpGs in classC",
          xlab = "Threshold (%)"
          #legend.title = "MPG Group"
) + ggplot2::geom_text(data=summary.table, aes(x = Threshold, 
                                               y = No_CpGs_classC,
                                               label=No_CpGs_classC), vjust=0)
dev.off()



pdf(paste0(class.dir, "classD.statistics.barplot", Sys.Date(), ".pdf"), height=4, width=5)
ggbarplot(summary.table, x = "Threshold", 
          y = "No_CpGs_classD",
          fill = "brown1",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          #palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          #sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "No. CpGs in classD",
          xlab = "Threshold (%)"
          #legend.title = "MPG Group"
) + ggplot2::geom_text(data=summary.table, aes(x = Threshold, 
                                               y = No_CpGs_classD,
                                               label=No_CpGs_classD), vjust=0)


dev.off()

