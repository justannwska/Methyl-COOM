## 04_Determine ClassA,B.R
## 
## These events are calculated on B cell specific CpGs
##
## R version 3.5.0
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
relative.methylation.changes <- readRDS(paste0(data.dir, "/Relative.meth.linear.model.on.B.cell.specific.CpGs.CLLs-2020-02-01.RDS"))





#######    Determine aberrant disease-specific methylation events (with a cutoff of 20%)  ####################################
# frequencies of methylation change more than 20%  methylation events (for classB)
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



# optional: plotting a histogram for methylation changes
# plot(density(frequencies.true$nr.patients))
# hist(frequencies.true$nr.patients, xlab= "Number of patients", ylab="Number of CpG sites", main="", 
#      cex=1, cex.lab=1, ylim=c(0,10000), las=2)
# abline(v=31, col="blue", lwd=2)

sample.fraction <- quantile(frequencies.true$nr.patients)





#  frequencies methylation change less than 20% (for classA)
frequencies.true.less <- data.frame("nr.patients" = as.numeric())
for(element in 1:nrow(relative.methylation.changes)){
  #print(element)
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







# optional: plotting a histogram for methylation changes
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
write.csv(frequencies.true, paste0(class.dir, "/frequencies.differentiation.related.more.than.0.2 cutoff.aberrant.methylation.events.per.patient", Sys.Date(), ".csv"))
write.csv(frequencies.true.less, paste0(class.dir, "/frequencies.differentiation.related.less.than.0.2 cutoff.aberrant.methylation.events.per.patient", Sys.Date(), ".csv"))






####################################################################################
# select the CpG sites with 20%  cutoff for variable number of patients (90% - 10% of patients)
# Later on one can decide which cutoff for the number of patients to use
number.of.patients <- ncol(relative.methylation.changes)

## determine an absolute number of patients belonging to each threshold
patient.fraction <- c( 0.9, 0.8, 0.75, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1)
#patient.fraction <- c(0.75)

patient.subset <- round(number.of.patients * patient.fraction)
  

## looping over all frequencies
classA.list <- list()
classB.list <- list()

relative.methylation.classA <- list()
relative.methylation.classB <- list()


#for(element in 1:6) {
for(element in 1:length(patient.subset)) {
  
  element.no <- patient.subset[[element]]
  
  ## determine number of CpGs in classA & B
  classA.list[[element]] <- data.frame(frequencies.true.less[frequencies.true.less$nr.patients >= as.numeric(element.no), ])
  classA.list[[element]]$cg <- as.character(classA.list[[element]]$cg)
  
  classB.list[[element]] <- data.frame(frequencies.true[frequencies.true$nr.patients >= as.numeric(element.no), ])
  classB.list[[element]]$cg <- as.character(classB.list[[element]]$cg)
  
  
  
  ## extract relative methylation values for the classA & B events
  
  # classA
  relmeth.objectA <- relative.methylation.changes[rownames(relative.methylation.changes) %in% classA.list[[element]]$cg, ]
  relmeth.objectA$cpg <- rownames(relmeth.objectA)
  
      if( all(relmeth.objectA$cpg %in% classA.list[[element]]$cg) ) {
        print("")
      } else {
        print("error: relative methylation was not extracted for all CpGs from classA")
      }
  
  # barplot per patient all CpGs
  data.for.barplot <- melt(relmeth.objectA[1:4,], id="cpg")
  
  p = ggplot(data = data.for.barplot, aes(x = cpg, y = value))
  p = p + geom_bar(stat='identity', width = 0.1) + geom_hline(yintercept= c(-20, 20), color="red")
  p = p + facet_grid(~variable)   
  p
  
  ggsave(paste0(class.dir, "//classA.perCpG.barplot", patient.fraction[element], Sys.Date(), ".pdf"),  plot= p, width= 30)
  
  
  
  # classB
  relmeth.objectB <- relative.methylation.changes[rownames(relative.methylation.changes) %in% classB.list[[element]]$cg, ]
  relmeth.objectB$cpg <- rownames(relmeth.objectB)
  
      if( all(relmeth.objectB$cpg %in% classB.list[[element]]$cg) ) {
        print("")
      } else {
        print("error: relative methylation was not extracted for all CpGs from classB")
      }
  
  # barplot per patient all CpGs
  data.for.barplot <- melt(relmeth.objectB[1:4,], id="cpg")
  
  p = ggplot(data = data.for.barplot, aes(x = cpg, y = value))
  p = p + geom_bar(stat='identity', width = 0.1) + geom_hline(yintercept= c(-20, 20), color="red")
  p = p + facet_grid(~variable)   
  p
  
  ggsave(paste0(class.dir, "//classB.perCpG.barplot", patient.fraction[element], Sys.Date(), ".pdf"),  plot= p, width= 30)
  
  
  # assign relative methylation to the list
  relative.methylation.classA[[element]] <- relmeth.objectA
  relative.methylation.classB[[element]] <- relmeth.objectB
  
  
  #check for duplicated CpGs (if duplicated an error will occur)
  #classA
    if(any(duplicated( relative.methylation.classA[[element]]$cpg)) )  {
      print("Error: Duplicated CpGs in extracted relative methylaion classA")
    } else {
      print("")
    }
  
  #classB
    if(any(duplicated(relative.methylation.classB[[element]]$cpg)) ) {
      print("Error: Duplicated CpGs in extracted relative methylaion classB")
    } else {
      print("")
    }
    
  
}

names(relative.methylation.classA) <- c( 0.9, 0.8, 0.75, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1)
names(relative.methylation.classB) <- c( 0.9, 0.8, 0.75, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1)
  

## saving a file with relative methylation for CpGs defined in classA and B
saveRDS(relative.methylation.classA, paste0(class.dir, "/classA_relative_methylation-", Sys.Date(), ".RDS"))
saveRDS(relative.methylation.classB, paste0(class.dir, "/classB_relative_methylation-", Sys.Date(), ".RDS"))





## extracting number of CLL-specific CpGs, depending on the threshold used 
CpG.count.classA <- sapply(relative.methylation.classA, function(get.cpg.count) { nrow(get.cpg.count) } )
CpG.count.classB <- sapply(relative.methylation.classB, function(get.cpg.count) { nrow(get.cpg.count) } )




## Creating a summary table with a number of events per class per number of patients
summary.table <- data.frame(patient.fraction*100, patient.subset, CpG.count.classA, CpG.count.classB)
colnames(summary.table) <- c("Threshold", "No_Patients", "No_CpGs_ClassA", "No_CpGs_ClassB")
summary.table[, c("No_CpGs_ClassA", "No_CpGs_ClassB")] <- apply(summary.table[, c("No_CpGs_ClassA", "No_CpGs_ClassB")], 2, function(x) as.numeric(as.character(x)));
summary.table$Threshold <- as.character(summary.table$Threshold)
write.xlsx(summary.table, paste0(class.dir, "/classA&B_statistics", Sys.Date(), ".xlsx"))




## plot the summary results
summary.data <- data.frame(rep(patient.fraction*100,2), c(CpG.count.classA, CpG.count.classB))
summary.data$class <- c( rep("classA", length(patient.subset)), rep("classB", length(patient.subset)) )
colnames(summary.data) <- c("Threshold", "No_CpGs", "Class")
summary.data$Threshold <- as.numeric(summary.data$Threshold)



# dotchart
ggdotchart(summary.data, x = "Threshold" , y = "No_CpGs",
           color = "Class",                                # Color by groups
           #palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
           #sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           group = "Class",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(summary.data$No_CpGs),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
)


# barplot
pdf(paste0(class.dir,"ClassA.statistics.barplot", Sys.Date(), ".pdf"), height=4, width=5)
ggbarplot(summary.table, x = "Threshold", 
          y = "No_CpGs_ClassA",
          fill = "dodgerblue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          #palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          #sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "No. CpGs in classA",
          xlab = "Threshold (%)"
          #legend.title = "MPG Group"
) + ggplot2::geom_text(data=summary.table, aes(x = Threshold, 
                                               y = No_CpGs_ClassA,
                                               label=No_CpGs_ClassA), vjust=0)
dev.off()



pdf(paste0(class.dir,"ClassB.statistics.barplot", Sys.Date(), ".pdf"), height=4, width=5)
ggbarplot(summary.table, x = "Threshold", 
          y = "No_CpGs_ClassB",
          fill = "brown1",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          #palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          #sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "No. CpGs in classB",
          xlab = "Threshold (%)"
          #legend.title = "MPG Group"
) + ggplot2::geom_text(data=summary.table, aes(x = Threshold, 
                                          y = No_CpGs_ClassB,
                                          label=No_CpGs_ClassB), vjust=0)


dev.off()



