#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(clusterProfiler))

gene_id_to_kegg_identifier <- function(orgdb, gene_id, keyType) {
    if(keyType == 'GID') {
        kegg_identifier <- 'GID'
    } else {
        kegg_identifier <- 'ENTREZID'
    }
    mapping <- as.data.table(select(get(orgdb), keys=gene_id, columns=c(kegg_identifier, keyType), keytype=keyType))
    setnames(mapping, names(mapping), c('key', 'id'))
    mapping <- mapping[match(gene_id, key)]
    return(mapping)
}

parser <- ArgumentParser(description='Run clusterprofiler')
parser$add_argument('--dge', help='File of differential gene expression table [required]', required=TRUE)
parser$add_argument('--orgdb', help='Organism database e.g. org.Hs.eg.db [required]', required=TRUE)
parser$add_argument('--dge-gene-id-column', help='Column of gene_id in DGE table [%(default)s]', default='gene_id')
parser$add_argument('--dge-score-column', help='Rank genes according to this column[%(default)s]', default='log2FoldChange')
parser$add_argument('--dge-group-column', help='The DGE table has results for several contrasts concatenated and indexed by this column. Apply clusterProfiler to each contrast group. If NULL, the DGE table has only one contrast group [%(default)s]', default=NULL)
parser$add_argument('--keyType', help='Type of gene identifier linking dge table to the annotation database. Popular options: ENSEMBL, ENTREZID, GID [%(default)s]', default='ENSEMBL')
parser$add_argument('--gsea-output-tsv', help='File for gsea output table', required=TRUE)
parser$add_argument('--kegg-organism', help='Organism abbreviation for gseKEGG. See http://www.genome.jp/kegg/catalog/org_list.html [%(default)s]', default='hsa')
parser$add_argument('--kegg-output-tsv', help='File for kegg output table', required=FALSE)

parser$add_argument('--version', '-v', action= 'version', version='0.2.0')

xargs <- parser$parse_args()

suppressPackageStartupMessages(library(xargs$orgdb, character.only=TRUE))

dge <- fread(xargs$dge)

stopifnot('tmp_dge_group' %in% names(dge) == FALSE)
stopifnot(xargs$dge_gene_id_column %in% names(dge))
stopifnot(xargs$dge_score_column %in% names(dge))

if(is.null(xargs$dge_group_column)) {
    dge[, tmp_dge_group := 'tmp']
} else {
    stopifnot(xargs$dge_group_column %in% names(dge))
    setnames(dge, xargs$dge_group_column, 'tmp_dge_group')
}

gseaOut <- list()
keggOut <- list()
for(cntr in unique(dge$tmp_dge_group)) {
    full <- dge[tmp_dge_group == cntr, c(xargs$dge_gene_id_column, xargs$dge_score_column), with=FALSE]
    setnames(full, names(full), c('gene_id', 'score'))

    if(length(full$gene_id) != length(unique(full$gene_id))) {
        stop('Duplicate gene ids found')
    }
    geneList <- full$score
    names(geneList) <- full$gene_id
    geneList <- na.omit(geneList)
    geneList <- sort(geneList, decreasing = TRUE)
    
    # GSEA GO
    gsea <- gseGO(geneList=geneList, 
                 ont="ALL", 
                 keyType=xargs$keyType, 
                 minGSSize=10, 
                 maxGSSize=1000, 
                 pvalueCutoff=1, 
                 verbose=FALSE, 
                 OrgDb=get(xargs$orgdb), # use get because you pass the object orgdb not the string 'org.XYZ.eg.db'
                 pAdjustMethod="BH",
                 seed=1234)

    gsea <- as.data.table(gsea@result)
    gsea[, tmp_dge_group := cntr]
    gseaOut[[length(gseaOut) + 1]] <- gsea

    # GSEA KEGG
    # We need to convert user identifier to the identifier used by KEGG
    if(!is.null(xargs$kegg_output_tsv)) {
        geneList <- data.table(gene_id=names(geneList), score=geneList)
        mapping <- gene_id_to_kegg_identifier(xargs$orgdb, geneList$gene_id, xargs$keyType)
        geneList <- merge(geneList, mapping, by.x='gene_id', by.y='key', sort=FALSE)
        geneList <- setNames(geneList$score, geneList$id)

        kegg <- gseKEGG(geneList, xargs$kegg_organism, pvalueCutoff=1, minGSSize=10, maxGSSize=1000, seed=1234)

        kegg <- as.data.table(kegg@result)
        kegg[, tmp_dge_group := cntr]
        keggOut[[length(keggOut) + 1]] <- kegg
    }
}
gsea <- rbindlist(gseaOut)
if(is.null(xargs$dge_group_column)) {
    gsea[, tmp_dge_group := NULL]
} else {
    setnames(gsea, 'tmp_dge_group', xargs$dge_group)
}
gsea[, qvalues := NULL]
gsea[, leading_edge := NULL]
gsea[, core_enrichment := NULL]
fwrite(gsea, xargs$gsea_output_tsv, sep='\t')


if(!is.null(xargs$kegg_output_tsv)) {
    kegg <- rbindlist(keggOut)
    if(is.null(xargs$dge_group_column)) {
        kegg[, tmp_dge_group := NULL]
    } else {
        setnames(kegg, 'tmp_dge_group', xargs$dge_group)
    }
    kegg[, qvalues := NULL]
    kegg[, leading_edge := NULL]
    kegg[, core_enrichment := NULL]
    fwrite(kegg, xargs$kegg_output_tsv, sep='\t')
}
