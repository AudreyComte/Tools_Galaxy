#source("/media/audrey/data/galaxy/tools/newTools/Deseq2_Report/diffexp_reports_compilation_galaxy.R")
#source("/media/audrey/data/galaxy/tools/newTools/Deseq2_Report/data_quality_galaxy.R")


## @knitr diffexp
library(stringr)

filename <- RDS

lfcshrink_type <- param3
FCcutoff <- param4
pCutoff <- param5

## Ploduce the MA plot and volcanoplot with matrices

cat("We use ",lfcshrink_type ," lfcshrink type as adaptive t prior shrinkage estimator for ranking and visualization.\n")
cat("The Sample of reference for the pairwise comparison is ",ref_level,".\n")

# Automation of the comparisons with the reference condition established in (config.yaml for deseq2.R) 
comparison <- resultsNames(dds)
comparison_df <- as.data.frame(comparison)[-1,]
outputname <- as.character(sub("condition_", "", comparison_df[[1]]))
condition <- strsplit(comparison_df[[1]], split = "_vs_")[[1]]
condition <- as.character(sub("condition_", "", condition))

# Obtenir les résultats d'expression différentielle pour le contraste courant
results <- results(dds, contrast = c("condition", condition))

# Réduction de la variance pour la visualisation
rlog_results <- lfcShrink(dds, coef = comparison_df[[1]], res = results, type = lfcshrink_type)
rlog_results <- rlog_results[order(rlog_results$padj),]

name <- paste0(mutant_level, "_vs_", ref_level, "_all_gene_stats.tsv")
write(name, file = output_names_stats, append = TRUE)
output_file2 <- file.path(absolute_path2, name)
# write.table(as.data.frame(rlog_results), file = paste("/home/audrey/Documents/Diff_expression/", outputname, "_all_genes_stats.tsv", sep = ""), quote = FALSE, sep = "\t")
write.table(as.data.frame(rlog_results), file = output_file2, quote = FALSE, sep = "\t")

# Génération des graphiques

## MA-Plot
xlim <- c(500, 5000)
ylim <- c(-2, 2)
plotMA(rlog_results, xlim = xlim, ylim = ylim, main = comparison_df[[1]])
idx <- identify(rlog_results$baseMean, rlog_results$log2FoldChange)

## Volcano Plot
FCcutoff <- as.numeric(FCcutoff)
FC <- log2(FCcutoff)
pCutoff <- as.numeric(pCutoff)
p <- pCutoff

# Initialisation des couleurs des points pour le volcanoplot
keyvals <- rep('grey75', nrow(rlog_results))
names(keyvals) <- rep('NS', nrow(rlog_results))

keyvals[which(abs(rlog_results$log2FoldChange) > FC & rlog_results$padj > p)] <- 'grey50'
names(keyvals)[which(abs(rlog_results$log2FoldChange) > FC & rlog_results$padj > p)] <- 'log2FoldChange'

keyvals[which(abs(rlog_results$log2FoldChange) < FC & rlog_results$padj < p)] <- 'grey25'
names(keyvals)[which(abs(rlog_results$log2FoldChange) < FC & rlog_results$padj < p)] <- '-Log10Q'

keyvals[which(rlog_results$log2FoldChange < -FC & rlog_results$padj < p)] <- 'blue2'
names(keyvals)[which(rlog_results$log2FoldChange < -FC & rlog_results$padj < p)] <- 'Signif. down-regulated'
sdr <- subset(rlog_results, (rlog_results$log2FoldChange < -FC) & (rlog_results$padj < p))
#write.table(as.data.frame(sdr), file = paste("/home/audrey/Documents/Diff_expression/", outputname, "_signif-down-regulated.txt", sep = ""), quote = FALSE, sep = "\t")
name_down <- paste0(mutant_level, "_vs_", ref_level, "_signif_down_regulated.txt")
output_file3 <- file.path(absolute_path3, name_down)

write(name_down, file = output_names_file_down, append = TRUE)

write.table(as.data.frame(sdr), file = output_file3, quote = FALSE, sep = "\t")

keyvals[which(rlog_results$log2FoldChange > FC & rlog_results$padj < p)] <- 'red2'
names(keyvals)[which(rlog_results$log2FoldChange > FC & rlog_results$padj < p)] <- 'Signif. up-regulated'
sur <- subset(rlog_results, (rlog_results$log2FoldChange > FC) & (rlog_results$padj < p))
#write.table(as.data.frame(sur), file = paste("/home/audrey/Documents/Diff_expression/", outputname, "_signif-up-regulated.txt", sep = ""), quote = FALSE, sep = "\t")
name_up <- paste0(mutant_level, "_vs_", ref_level, "_signif_up_regulated.txt")
output_file4 <- file.path(absolute_path4, name_up)

write(name_up, file = output_names_file_up, append = TRUE)

write.table(as.data.frame(sur), file = output_file4, quote = FALSE, sep = "\t")


# Création du volcanoplot
volcanoplot_padj <- EnhancedVolcano(rlog_results,
lab = rownames(rlog_results),
x = 'log2FoldChange',
y = "padj",
xlim = c(-12, 12),
xlab = bquote(~Log[2]~ 'fold change'),
ylab = bquote(~-Log[10]~adjusted~italic(P)),
title = "",
subtitle = "",
axisLabSize = 12,
titleLabSize = 15,
pCutoff = p,
FCcutoff = FC,
pointSize = 1.5,
labSize = 2,
colCustom = keyvals,
colAlpha = 1,
legendPosition = 'bottom',
legendLabSize = 10,
legendIconSize = 1.5,
gridlines.major = FALSE,
gridlines.minor = FALSE,
border = 'full',
borderWidth = 0.7,
borderColour = 'black')
plot(volcanoplot_padj)


if(length(gene_name_list)!=0){
	for (i in 1:length(gene_name_list)) {
        # adjusted p-value can be NA (>1) so condition doesn't work
        gene_name=gene_name_list[[i]]
        foldc = rlog_results[which(row.names(rlog_results) == gene_name),"log2FoldChange"]
        pval = rlog_results[which(row.names(rlog_results) == gene_name),"padj"]
        if (!is.na(pval)){
            if ((foldc < -FC) && (pval < p)) {
                cat(paste("\n",gene_name," is differentially expressed (significantly down regulated).\n",sep=""))
            }
            if ((foldc > FC) && (pval < p)) {
                cat(paste("\n",gene_name," is differentially expressed (significantly up regulated).\n",sep=""))
            }
        }
    }
}


## @knitr toppadj

#stat_file=paste(mutant_level,"_vs_",ref_level,"_all_genes_stats.tsv", sep="")
nbpval=param7
mutant_level_file <- read_delim(output_file2, "\t", escape_double = FALSE, trim_ws = TRUE)

# Top adjusted p-value
top_adjpval <- mutant_level_file[1:nbpval,1]
# Keep the normalized count values of these genes
dds <- readRDS(RDS)
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))

normalized_counts <- rownames_to_column(normalized_counts, "Gene")
topnorm <- data.frame()
for (i in 1:nrow(top_adjpval)){
    genename=as.character(top_adjpval[i,1])
    ligne=subset(normalized_counts,normalized_counts$Gene==genename)
    topnorm <- rbind(topnorm,ligne)
}
topnorm_length=length(topnorm)
gathered_topnorm <- topnorm %>% gather(colnames(topnorm)[2:topnorm_length], key = "samplename", value = "normalized_counts")

# Counts colored by sample group
coldata_read <- read.delim(param8, header=TRUE, comment.char="#", quote="")
coldata <- coldata_read[,-3]
names(coldata)[names(coldata) == "project"] <- "samplename"
gathered_topnorm <- inner_join(coldata, gathered_topnorm)

cat(paste("\nWe plot the normalized count values for the top ",nbpval," differentially expressed genes (by padj values).\n",sep=""))

cat(paste("\nThese genes are selected with the padj-value obtened from the comparison ",mutant_level," vs ", ref_level,"\n", sep=""))

p=ggplot(gathered_topnorm) +
        geom_point(aes(x = Gene, y = normalized_counts, color = condition)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("Normalized Counts") +
        ggtitle(paste("Top ", nbpval, " Significant DE Genes",sep="")) +
        theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(plot.title = element_text(hjust = 0.5))
ggplotly(p)