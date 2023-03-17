#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description='Run GSEA on differential genes') 
parser$add_argument('--dge', help='Table of DGE [%(default)s]', default='-')
parser$add_argument('--geneid-column', help='Name of gene_id column in dge table [%(default)s]', default='gene_id')
parser$add_argument('--score-column', help='Name of score column in dge table [%(default)s]', default='log2FoldChange')
parser$add_argument('--ensembl-version', type='integer', help='Version of ensembl database to use [%(default)s]', default=NULL)
parser$add_argument('--biomart-dataset', help='Use this biomart dataset [%(default)s]', default='hsapiens_gene_ensembl')
parser$add_argument('--nproc', type='integer', help='Number of processors for fgsea [%(default)s]', default=16)
parser$add_argument('--out-tab', help='File for output table [%(default)s]', default='-')
xargs <- parser$parse_args()

# xargs <- list(dge='Homo_sapiens.GRCh38.105/RNAseq/deseq/dge.tsv.gz', biomart_dataset='hsapiens_gene_ensembl', ensembl_version=105, geneid_column='gene_id', score_column='log2FoldChange')

list_ensembl_go <- function(ensembl_mart, genes) {
    genes <- unique(genes)
    chunks <- split(genes, ceiling(seq_along(genes)/1000))
    go_terms <- list()
    i <- 1
    for(genes in chunks) {
        write(sprintf('Retrieving GO terms | chunk %s/%s', i, length(chunks)), stderr())
        go_terms[[length(go_terms) + 1]] <- getBM(attributes=c('ensembl_gene_id', 'go_id'), filters=c('ensembl_gene_id'), values=list(genes), mart=ensembl_mart)
        i <- i+1
    }
    go_terms <- rbindlist(go_terms)
    go_terms <- go_terms[go_id != '']

    gene2go <- tapply(go_terms$go_id, go_terms$ensembl_gene_id, FUN=c)
    return(gene2go)
}

full_go_to_gene <- function(gene2go) {
    # Prepare list of GO terms and assiociated genes suitable for fgsea
    #
    # gene2go: list of gene_id (key) and associated vector of go terms (value). 
    # Return: list of GO terms (key) and associated vector of genes. The vector
    # of genes is augmented to include all the ancestors terms.
    goa <- list()
    for(go in c(GOMFANCESTOR, GOBPANCESTOR, GOCCANCESTOR)) {
        xx <- as.list(go)
        goa <- c(goa, xx)
    }

    gene2fullgo <- list()
    for(gene in names(gene2go)) {
        go <- gene2go[[gene]]
        ancestors <- unlist(goa[go])
        fullgo <- unique(c(go, ancestors))
        fullgo <- fullgo[fullgo != 'all']
        gene2fullgo[[gene]] <- fullgo
    }
    gene2fullgo <- data.table(ensembl_gene_id=rep(names(gene2fullgo), sapply(gene2fullgo, length)), go_id=unlist(gene2fullgo))
    go2gene <- tapply(gene2fullgo$ensembl_gene_id, gene2fullgo$go_id, FUN=c)
    return(go2gene)

    # Credit: https://stackoverflow.com/questions/73634484/efficiently-convert-two-columns-data-table-to-list/73634548#73634604
    # go_list <- tapply(go_terms$ensembl_gene_id, go_terms$go_id, FUN=c)
    # return(go_list)
}

if(xargs$dge == '-') {
    fin <- 'file:///dev/stdin'
} else {
    fin <- xargs$dge
}
dge <- fread(fin, select=c(xargs$geneid_column, xargs$score_column), col.names=c('gene_id', 'log2FoldChange'))
dge <- unique(dge)
stopifnot(length(dge$gene_id) == length(unique(dge$gene_id)))

# Get GO annotation 
ensembl <- useEnsembl(biomart = 'genes',
                      dataset = xargs$biomart_dataset,
                      version = xargs$ensembl_version)

scores <- dge$log2FoldChange
names(scores) <- dge$gene_id
scores <- sort(scores, decreasing=TRUE)
gene2go <- list_ensembl_go(ensembl, genes=names(scores))
go2gene <- full_go_to_gene(gene2go)

write('Running fgsea...', stderr())
gsea <- fgseaMultilevel(go2gene, stats=scores, minSize=10, scoreType='std', nproc=xargs$nproc)
setnames(gsea, 'pathway', 'go_id')

goterms <- Term(GOTERM)
goterms <- data.table(go_id=names(goterms), go_name=goterms)

gsea <- merge(gsea, goterms, by='go_id')
# idx <- which(names(gsea) != 'leadingEdge') # Move leadingEdge col to last position
# setcolorder(gsea, idx)
gsea[, leadingEdge := NULL]
gsea <- gsea[order(-abs(NES))]

if(xargs$out_tab == '-') {
    fout <- ''
} else {
    fout <- xargs$out_tab
}
fwrite(gsea, fout, sep='\t', na='NA', quote=FALSE)
