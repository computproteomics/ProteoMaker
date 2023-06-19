### needs allBs and BenchMatrix 
library(reshape2)
library(ggplot2)
library(lattice)

# take specific parameter for proportion of miscleavages
prop_mcs <- 0.1

dat_with_prop_mcs <- rownames(BenchMatrix)[which(BenchMatrix$PropMissedCleavages == prop_mcs)[1]]

filename <- paste0(resultFilePath,"/outputDataAnalysis_",dat_with_prop_mcs,".RData")
load(filename)

StatsPep$MC <- sapply(StatsPep$MC, unique)
melted <- melt(StatsPep[,c("MC",Param$QuantColnames)], id.vars="MC")

ggplot(melted, aes(x=value, fill=MC)) +
  geom_histogram(binwidth=0.25, aes(x=value,y=..density..)) + 
  facet_grid(MC~.) + xlab("abundance") + title(paste("Proportion of",prop_mcs, "miscleavages"))


barplot(t(BenchMatrix[1:10,paste0("propMisCleavedPeps.",0:6)]), beside = T, names.arg = BenchMatrix$PropMissedCleavages[1:10],
        xlab="Proportion miscleavages")

barplot(t(BenchMatrix[1:10,"propMisCleavedProts"]), beside = T, names.arg = BenchMatrix$PropMissedCleavages[1:10],
        xlab="Proportion miscleavages", main="proteins with miscleaved peptides")


## checking whether quantNoise and PropMissedCleavages have some combined effect
mc_map <- matrix(NA, nrow=length(unique(unique(BenchMatrix$PropMissedCleavages))),  
                 ncol=length(unique(BenchMatrix$QuantNoise)), dimnames=list(x=unique(BenchMatrix$PropMissedCleavages), y=unique(BenchMatrix$QuantNoise)))
for (i in 1:nrow(BenchMatrix))
  mc_map[BenchMatrix$PropMissedCleavages[i], BenchMatrix$QuantNoise[i]] <- BenchMatrix$numQuantProtGroups[i]


levelplot(mc_map, xlab="PropMissedCleavages", ylab="QuantNoise")





par(mfrow=c(1,1))
hist(Stats$C_1_R_5, 100)
