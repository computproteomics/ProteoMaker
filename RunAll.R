################################################################################
#                     TO TUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

#####################
## Load parameters
#####################
source("Parameter.R")
#####################

#####################
## Run the sample preparation simulation:
#####################
source("01_GenerateGroundTruth.R")
# Create the initial list of proteoforms:
proteoforms <- samplePreparation(fasta.path = Param$PathToFasta, parameters = Param)
# Create the full structure of proteoforms along with abundances and ground truth expression patterns:
GroundTruth <- addProteoformAbundance(proteoforms = proteoforms, parameters = Param)
rm(proteoforms)
# # Save GroundTruth for analysis:
# save(GroundTruth, file = "RData/GroundTruthAbs.RData")
#####################

# #####################
# ## Visualisation Ground truth data:
# #####################
# library(ggplot2)
# library(reshape2)
# # Identifications
# cat("Number of unique protein accessions:", length(unique(GroundTruth$Accession)), "\n")
# gtab <- GroundTruth
# gtab$Regulated <- !is.na(gtab$Regulation_Amplitude)
# if (Param$FracModProt > 0) {
#   gtab$modProt <- !(sapply(gtab$PTMType, is.null))
#   cat("Number of proteoforms per accession:\n")
#   hist(table(GroundTruth$Accession), col = "cornflowerblue", xlab = "Number of proteoforms per accession", main = "", breaks = 20)
#   cat("Number of protein accessions with only one proteoform:", sum(table(GroundTruth$Accession) == 1), "\n")
#   cat("Parameter to select the fraction of proteins to be modified:", Param$FracModProt)
#   cat("Number of accessions with modifications:", length(unique(gtab$Accession[gtab$modProt])))
#   cat("Number of regulated proteoforms:", sum(gtab$Regulated))
#   cat("Number of regulated accessions:", length(unique(gtab$Accession[gtab$Regulated])))
#   
#   gtab2 <- melt(table(gtab$Accession, gtab$modProt))
#   names(gtab2)[2] <- "isModified"
# 
#   g <- ggplot(data = gtab2, aes(x = value, fill = isModified)) +
#     geom_bar(position = "dodge", col = "black") +
#     scale_fill_manual(values = c("grey30", "#DF8F44FF")) +
#     theme_bw() +
#     xlab("Number of proteoforms per protein")
#   print(g)
#   venn::venn(list("modified" = unique(gtab$Accession[gtab$modProt]), "unmodified" = unique(gtab$Accession[!(gtab$modProt)])))
#   
#   gtab2 <- melt(table(gtab$Accession[gtab$modProt], gtab$Regulated[gtab$modProt]))
#   names(gtab2)[2] <- "isRegulated"
#   g <- ggplot(data = gtab2, aes(x = value, fill = isRegulated)) +
#     geom_bar(position = "dodge", col = "black") +
#     scale_fill_manual(values = c("grey30", "#c10534")) +
#     theme_bw() +
#     xlab("Number of proteoforms per protein") +
#     labs(title = "Only modified accessions")
#   print(g)
#   gtab2 <- melt(table(gtab$Accession[!(gtab$modProt)], gtab$Regulated[!(gtab$modProt)]))
#   names(gtab2)[2] <- "isRegulated"
#   g <- ggplot(data = gtab2, aes(x = value, fill = isRegulated)) +
#     geom_bar(position = "dodge", col = "black") +
#     scale_fill_manual(values = c("grey30", "#c10534")) +
#     theme_bw() +
#     xlab("Number of proteoforms per protein") +
#     labs(title = "Only unmodified accessions")
#   print(g)
#   
#   # Number of sites
#   par(mar = c(6,3,1,1))
#   modtab <- gtab[gtab$modProt,]
#   modSites <- sapply(seq_len(nrow(modtab)), function(x) {
#     paste(modtab$Accession[x], modtab$PTMType[[x]], modtab$PTMPos[[x]])
#   })
#   modSites <- unlist(modSites)
#   cat("There are", length(unique(modSites)), "unique modifications identified (accession + localised modif.).\n")
#   hist(table(modSites), col = "cornflowerblue", xlab = "Number occurence of a given modification (accession + localised modif.)", main = "", breaks = 20)
# }
# 
# # Quantification
# par(mar = c(6,3,1,1))
# hist(unlist(sapply(seq_along(gtab$Regulation_Amplitude), function(x) gtab$Regulation_Amplitude[x]*gtab$Regulation_Pattern[[x]])), 
#      col = "#c10534", 
#      main = "", 
#      xlab = "Amplitude of regulation", 
#      breaks = 50)
# cat("There are", sum(is.na(gtab[,grepl("^C_", names(gtab))])), "missing values in the table.\n")
# 
# if (Param$FracModProt > 0) {
#   gtab2 <- melt(gtab[,grepl("^C_|Regulated|modProt", names(gtab))], id.vars = c("Regulated", "modProt"))
#   g <- ggplot(gtab2, aes(x = variable, y = value, fill = Regulated)) +
#     geom_violin(alpha = 0.3) +
#     scale_fill_manual(values = c("grey30", "#c10534")) +
#     geom_boxplot(width = 0.3, position = position_dodge(width = 1)) +
#     theme_bw()
#   print(g)
#   names(gtab2)[names(gtab2) == "modProt"] <- "isModified"
#   g <- ggplot(gtab2, aes(x = variable, y = value, fill = isModified)) +
#     geom_violin(alpha = 0.3) +
#     scale_fill_manual(values = c("grey30", "#DF8F44FF")) +
#     geom_boxplot(width = 0.3, position = position_dodge(width = 1)) +
#     theme_bw()
#   print(g)
# } else {
#   gtab2 <- melt(gtab[,grepl("^C_|Regulated|modProt", names(gtab))], id.vars = c("Regulated"))
#   g <- ggplot(gtab2, aes(x = variable, y = value, fill = Regulated)) +
#     geom_violin(alpha = 0.3) +
#     scale_fill_manual(values = c("grey30", "#c10534")) +
#     geom_boxplot(width = 0.3, position = position_dodge(width = 1)) +
#     theme_bw()
#   print(g)
# }
# #####################

#####################
## Digestion and sample peptide enrichment:
#####################
source("02_Digestion.R")
# Digest all the proteoforms and get peptide table:
peptable <- DigestGroundTruth(GroundTruth = GroundTruth, parameters = Param)
peptable <- mapQuanToDigestionProd(DigestedProt = peptable)
BeforeMS <- filterDigestedProt(peptable, Param)
# # Save peptable before filter for analysis:
save(BeforeMS, file = "RData/BeforeMS.RData")
#####################