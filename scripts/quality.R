library(data.table)
library(ggplot2)
library(DESeq2)
library(ggrepel)

ss <- fread(snakemake@input[['ss']])
ss[, patient_id := as.factor(as.character(patient_id))]
ss[, neutrophils := as.factor(neutrophils)]
fin <- snakemake@input[['cnt']]  # Sys.glob('Homo_sapiens.GRCh38.105/RNAseq/featureCounts/*.mapq0.tsv')

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

fwrite(cnt, snakemake@output[['cnt']], sep= '\t')

biotypes <- cnt[, list(n_reads= sum(count), n_genes= .N, n_expr_gt10= sum(count > 10), n_expr_gt100= sum(count > 100)),
    by= list(gene_biotype, library_id)]

biotypes <- biotypes[, list(gene_biotype, n_reads, n_genes, n_expr_gt10, n_expr_gt100, pct_reads= 100 * n_reads/sum(.SD$n_reads)), 
    by= library_id]
biotypes <- biotypes[order(, library_id, -pct_reads)]
fwrite(biotypes, snakemake@output[['biotypes']], sep= '\t')

dcnt <- dcast(data= cnt[gene_biotype %in% c('protein_coding')], Geneid ~ library_id, 
    value.var= 'count', fill= 0)

dcnt <- as.matrix(dcnt, rownames= 'Geneid')
dcnt <- dcnt[, match(ss$library_id, colnames(dcnt))]
deseq <- DESeqDataSetFromMatrix(dcnt, colData= ss, design= ~patient_id + neutrophils)
vsd <- assay(vst(deseq, blind= TRUE))

lvsd <- melt(as.data.table(vsd, keep.rownames='gene_id'), id.vars='gene_id', variable.name='library_id', value.name='gex')
fwrite(lvsd, snakemake@output[['vsd']], sep='\t')

select <- order(apply(vsd, 1, var), decreasing=TRUE)[1:1000]
pca <- prcomp(t(vsd[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
d <- data.table(PC1= pca$x[,1], PC2= pca$x[,2], library_id= colnames(vsd))
pcaout <- merge(d, ss[, list(library_id, patient_id, neutrophils)], by= 'library_id', all= TRUE)
pcaout[, neutrophils := gsub('_', ' ', neutrophils)]
xy <- range(d[, list(PC1, PC2)])
gg1 <- ggplot(data= pcaout, aes(x= PC1, y= PC2, colour=neutrophils, label= patient_id)) +
    geom_point() +
    geom_text_repel(show.legend=FALSE) +
    scale_color_brewer(palette='Set1') +
    ylim(xy[1], xy[2]) +
    xlim(xy[1], xy[2]) +
    xlab(sprintf('PC1: %.1f%%', percentVar[1] * 100)) +
    ylab(sprintf('PC2: %.1f%%', percentVar[2] * 100)) +
    ggtitle('PCA of gene expression') +
    theme_light() 
ggsave(snakemake@output[['pca']], width= 14, height= 10, units= 'cm')

cormat <- as.data.table(cor(vsd), keep.rownames= 'library_1')
cormat <- melt(cormat, id.vars= 'library_1', variable.name= 'library_2', value.name= 'cor')
cormat[, library_2 := as.character(library_2)]
cormat[, library_2 := factor(library_2, sort(unique(library_2), decreasing= TRUE))]
gg <- ggplot(data= cormat, aes(x= library_1, y= library_2, fill= cor, label= sprintf('%.2f', cor))) +
    geom_tile() +
    geom_text(size= 3) +
    xlab('') +
    ylab('') +
    ggtitle('Correlation in gene expression')
ggsave(snakemake@output[['corr']], width= 14, height= 10, units= 'cm')

lvsd <- melt(as.data.table(vsd, keep.rownames= 'gene_id'), id.vars= 'gene_id', value.name= 'gex', variable.name= 'library_id')
lvsd <- merge(lvsd, ss[, list(library_id, RIN)], by= 'library_id')
lvsd[, library_rin := sprintf('%s | RIN: %s', library_id, RIN)]
xord <- unique(lvsd[order(RIN, library_id)]$library_rin)
lvsd[, library_rin := factor(library_rin, xord)]
gg <- ggplot(data= lvsd, aes(x= gex)) +
    geom_histogram(bins= 25) +
    facet_wrap(~library_rin) +
    xlab('Normalized expression') +
    ylab('Number of genes') +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'))
ggsave(snakemake@output[['gex']], width= 18, height= 16, units= 'cm')


smry <- list()
for(x in fin){
    dat <- fread(paste0(x, '.summary'))
    smry_col <- names(dat)[ncol(dat)]
    library_id <- sub('\\.bam$', '', basename(smry_col))
    setnames(dat, names(dat), c('category', 'count'))
    dat[, library_id := library_id]
    smry[[length(smry) + 1]] <- dat
}
smry <- rbindlist(smry)
smry <- smry[category != 'Unassigned_Secondary' & count > 0]
smry <- smry[, list(category, count, pct= 100 * count / sum(.SD$count)), library_id]
fwrite(smry, snakemake@output[['count_summary']], sep= '\t')
