#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggrepel))

biomart_gene_description <- function(genes, ensembl_version, biomart_dataset='hsapiens_gene_ensembl') {
    mart <- useEnsembl(biomart = 'genes',
                          dataset = biomart_dataset,
                          version = ensembl_version)
    genes <- data.table(gene_id=unique(genes))
    chunks <- split(genes$gene_id, ceiling(seq_along(genes$gene_id)/1000))
    biomart_out <- list()
    i <- 1
    for(g in chunks) {
        write(sprintf('Retrieving genes | chunk %s/%s', i, length(chunks)), stderr())
        biomart_out[[length(biomart_out) + 1]] <- getBM(attributes=c('ensembl_gene_id', 'description'), filters=c('ensembl_gene_id'), values=list(g), mart=mart)
        i <- i+1
    }
    biomart_out <- rbindlist(biomart_out)
    setnames(biomart_out, 'ensembl_gene_id', 'gene_id')
    biomart_out <- merge(genes, biomart_out, by='gene_id', all.x=TRUE)
    stopifnot(nrow(biomart_out) == nrow(genes))
    return(biomart_out)
}

parser <- ArgumentParser(description='Run DEseq2') 
parser$add_argument('--sample-sheet', required=TRUE)
parser$add_argument('--count-files', nargs='+', help='Count files from featureCounts; one for each library in sample sheet. library_id extracted from header', required=TRUE)
parser$add_argument('--ensembl-version', type='integer', help='Version of ensembl database to use [%(default)s]', default=NULL)
parser$add_argument('--biomart-dataset', help='Use this biomart dataset [%(default)s]', default='hsapiens_gene_ensembl')
parser$add_argument('--out-dge', help='Output file for DGE', required=TRUE)
parser$add_argument('--out-maplot', help='Output file for maplot', required=TRUE)
parser$add_argument('--out-boxplot', help='Output file for boxplot of log2 fold-changes', required=TRUE)

xargs <- parser$parse_args()

ss <- fread(xargs$sample_sheet)
ss[, patient_id := as.factor(as.character(patient_id))]
ss[, neutrophils := factor(neutrophils, c('normal_density', 'low_density'))]
fin <- xargs$count_files

cnt <- list()
for(x in fin){
    dat <- fread(x)
    cnt_col <- names(dat)[ncol(dat)]
    library_id <- sub('\\.bam$', '', basename(cnt_col))
    dat <- dat[, c('Geneid', 'Length', 'gene_name', 'gene_biotype', cnt_col), with= FALSE]
    dat[, library_id := library_id]
    setnames(dat, cnt_col, 'count')
    cnt[[length(cnt) + 1]] <- dat
}
cnt <- rbindlist(cnt)

dcnt <- dcast(data= cnt[gene_biotype %in% c('protein_coding', 'lncRNA')], Geneid ~ library_id, 
    value.var= 'count', fill= 0)

dcnt <- as.matrix(dcnt, rownames= 'Geneid')
dcnt <- dcnt[, match(ss$library_id, colnames(dcnt))]
ss[, cRIN := RIN - mean(ss$RIN)]
deseq <- DESeqDataSetFromMatrix(dcnt, colData= ss, design= ~patient_id + neutrophils)
deseq <- DESeq(deseq)

dge <- list()
for(cntr in c('neutrophils_low_density_vs_normal_density')){
    res <- lfcShrink(deseq, coef= cntr, type= 'ashr')
    res <- as.data.table(as.data.frame(res), keep.rownames= 'gene_id')
    res[, contrast := cntr]
    dge[[length(dge) + 1]] <- res
}
dge <- rbindlist(dge)
# dge <- dge[!is.na(padj)]

genes <- unique(cnt[, list(gene_id= Geneid, gene_name, gene_biotype)])
dge <- merge(dge, genes, by= 'gene_id')

desc <- biomart_gene_description(dge$gene_id, xargs$ensembl_version, xargs$biomart_dataset)
dge <- merge(dge, desc, by='gene_id', sort=FALSE)

fwrite(dge[order(padj)], xargs$out_dge, sep= '\t')

# MAPLOT

dge[, coldens := densCols(log10(baseMean), log2FoldChange, 
    nbin= 256, bandwidth= 2/3, colramp = colorRampPalette(rev(rainbow(10)))), by= contrast]
dge[, contrast_name := gsub('_', ' ', contrast)]

lfc <- 1
n_dge <- rbind(
    dge[, list(direction= 'Up', N= sum(padj < 0.01 & log2FoldChange > lfc)), by= contrast_name],
    dge[, list(direction= 'Down', N= sum(padj < 0.01 & log2FoldChange < -lfc)), by= contrast_name]
)

options(scipen= 9)
gg <- ggplot(data= dge, aes(baseMean, log2FoldChange, label= gene_name)) +
    geom_point(cex= 0.1, aes(colour= coldens)) +
    geom_text(data= n_dge[direction == 'Up'], x= -Inf, y= Inf, hjust= -0.1, vjust= 1.5, aes(label= sprintf('Up: %s', N))) +
    geom_text(data= n_dge[direction == 'Down'], x= -Inf, y= -Inf, hjust= -0.1, vjust= -0.5, aes(label= sprintf('Down: %s', N))) +
    scale_color_identity() +
    scale_x_log10() +
    scale_y_continuous(breaks= scales::pretty_breaks(n= 10)) +
    geom_hline(yintercept= 0, colour= 'black', linetype= 'dashed') +
    geom_hline(yintercept= c(lfc, -lfc), colour= 'grey30', linetype= 'dashed') +
    facet_wrap(~contrast_name, ncol= 2) +
    ylab('log2 fold-change') +
    xlab('Average expression') +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'))
ggsave(xargs$out_maplot, width= 12, height= 10, units= 'cm')

dge[, upregulated_in := ifelse(log2FoldChange > 0, 'Low density', 'High density')]
low <- dge[padj < 0.01 & log2FoldChange < 0][order(-abs(log2FoldChange))][1:10]
low[, gene_name := ifelse(gene_name == '', gene_biotype, gene_name)]
high <- dge[padj < 0.01 & log2FoldChange > 0][order(-abs(log2FoldChange))][1:10]
high[, gene_name := ifelse(gene_name == '', gene_biotype, gene_name)]

gg <- ggplot(data=dge[padj < 0.01], aes(x=upregulated_in, y=abs(log2FoldChange))) +
    geom_quasirandom(cex=0.1, colour='grey30') +
    geom_point(data=dge[padj < 0.01, list(log2FoldChange=median(abs(log2FoldChange))), by=upregulated_in], cex=10, pch='-', colour='red') +
    geom_text_repel(data=low, xlim=c(NA, 1.5), aes(label=gene_name), cex=2) +
    geom_text_repel(data=high, xlim=c(1.5, NA), aes(label=sub(' \\[.*', '', description)), cex=2) +
    ylab('log2 fold-change') +
    xlab('Genes upregulated in') +
    scale_y_continuous(breaks= scales::pretty_breaks(n= 10)) +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'))
ggsave(xargs$out_boxplot, width= 10, height= 10, units= 'cm')


## Effect of RIN (not used):
# deseq <- DESeqDataSetFromMatrix(dcnt, colData= ss, design= ~ cRIN)
# deseq <- DESeq(deseq)
# res <- lfcShrink(deseq, coef= 'cRIN', type= 'ashr')
# res <- as.data.table(as.data.frame(res), keep.rownames= 'gene_id')
# res[order(pvalue)]

