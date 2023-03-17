rule deseq:
    input:
        ss= config['sample_sheet'],
        cnt= expand('{species}/{library_type}/featureCounts/{library_id}.mapq10.tsv', zip, 
                library_type= ss[ss.library_type == 'RNAseq'].library_type, 
                species= ss[ss.library_type == 'RNAseq'].species, 
                library_id= ss[ss.library_type == 'RNAseq'].library_id),
    output:
        dge= '{species}/{library_type}/deseq/dge.tsv.gz',
        maplot= '{species}/{library_type}/deseq/maplot.pdf',
        boxplot= '{species}/{library_type}/deseq/boxplot.pdf',
    params:
        ensembl_version= lambda wc: species[species.species == wc.species].ensembl_version.iloc[0],
    shell:
        r"""
        {workflow.basedir}/scripts/deseq.r --sample-sheet {input.ss} --count-files {input.cnt} \
            --out-dge {output.dge} --out-maplot {output.maplot} --out-boxplot {output.boxplot} \
            --ensembl-version {params.ensembl_version}
        """


rule clusterprofiler:
    input:
        dge= '{species}/{library_type}/deseq/dge.tsv.gz',
    output:
        gsea= '{species}/{library_type}/clusterprofiler/gsea.tsv.gz',
        # Not working now due to https://github.com/YuLab-SMU/clusterProfiler/issues/561
        # kegg= '{species}/{library_type}/clusterprofiler/kegg.tsv.gz',
    shell:
        r"""
        {workflow.basedir}/scripts/clusterprofiler.r --dge {input.dge} \
                --orgdb org.Hs.eg.db \
                --dge-gene-id-column gene_id \
                --dge-score-column log2FoldChange \
                --dge-group-column contrast \
                --keyType ENSEMBL \
                --gsea-output-tsv {output.gsea}
        """


# rule gsea:
#     input:
#         dge= '{species}/{library_type}/deseq/dge.tsv.gz',
#     output:
#         gsea= '{species}/{library_type}/gsea/{contrast}.gsea.tsv.gz',
#     params:
#         ensembl_version= lambda wc: species[species.species == wc.species].ensembl_version.iloc[0],
#     shell:
#         r"""
#         zcat {input.dge} | grep -P -w 'gene_id|{wildcards.contrast}' \
#         | {workflow.basedir}/scripts/gsea.r --ensembl-version {params.ensembl_version} --out-tab {output.gsea}
#         """


rule parse_msigdb:
    input:
        yaml=os.path.join(workflow.basedir, 'data/hif1a.msigdb.yaml'),
    output:
        tsv=temp('hif1a.msigdb.tsv'),
    run:
        import yaml
        with open(input.yaml, 'r') as f:
            y = yaml.safe_load(f)
        gs = y['SEMENZA_HIF1_TARGETS']['geneSymbols']
        tsv = pandas.DataFrame({'source':'msigdb', 'gene_name':'HIF1A', 'targets':gs})
        tsv.to_csv(output.tsv, sep='\t', index=False)


rule target_profiles:
    input:
        dge='{species}/{library_type}/deseq/dge.tsv.gz',
        targets='hif1a.msigdb.tsv',
        gff='{species}/ref/annotation.gtf',
    output:
        rank='{species}/{library_type}/misc/hif1a_msigdb_targets_rank.pdf',
        bees='{species}/{library_type}/misc/hif1a_msigdb_targets_bees.pdf',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)

dge <- fread('{input.dge}')
targets <- fread('{input.targets}')
gff <- fread(cmd="awk '$3 == \"gene\"' {input.gff} | grep 'gene_name'", header=FALSE)

gff[, gene_id := sub('".*', '', sub('.*gene_id "', '', V9))]
gff[, gene_name := sub('".*', '', sub('.*gene_name "', '', V9))]
stopifnot(targets$targets %in% gff$gene_name)

targets <- merge(targets, gff[, list(gene_id, gene_name)], by.x='targets', by.y='gene_name')
dge[, is_target := (gene_id %in% targets$gene_id) | gene_name == 'HIF1A']
dge[, rank := rank(log2FoldChange), by=contrast]

gg <- ggplot(data=dge, aes(x=rank, y=log2FoldChange)) +
    geom_hline(yintercept=0, colour='black', linetype='dashed', size=0.1) +
    geom_point(cex=0.25, colour='grey60') +
    geom_point(data=dge[is_target == TRUE], cex=0.5, colour='red') +
    geom_text_repel(data=dge[is_target == TRUE], aes(label=gene_name), colour='black', size=2.5, max.overlaps=20, segment.color='grey80') +
    facet_wrap(~contrast) +
    xlab('Gene rank (ranked on fold-change)') +
    ylab('Log2 fold-change') +
    ggtitle('Differential expression of HIF1A and its targets from MSigDB') +
    theme_light() +
    theme(strip.text=element_text(colour='black'))
ggsave('{output.rank}', width=14, height=12, units='cm')

dge[, direction := ifelse(log2FoldChange > 0, sub('_vs_', ' > ', contrast), sub('_vs_', ' < ', contrast))]
dge[, direction := gsub('neutrophils_|_density', '', direction)]

gg <- ggplot(data=dge, aes(x=direction, y=abs(log2FoldChange))) +
    geom_quasirandom(colour='grey60', cex=0.25) +
    geom_point(data=dge[is_target == TRUE], cex=0.5, colour='red') +
    geom_text_repel(data=dge[is_target == TRUE], aes(label=gene_name), colour='black', size=2.5, max.overlaps=20, segment.color='grey80') +
    facet_wrap(~contrast) +
    xlab('') +
    ylab('Log2 fold-change') +
    ggtitle('Differential expression of HIF1A and its targets from MSigDB') +
    theme_light() +
    theme(strip.text=element_text(colour='black'), plot.title = element_text(size=10), axis.text.x=element_text(size=12, colour='black'))
ggsave('{output.bees}', width=12, height=12, units='cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """
        
rule volcano:
    input:
        dge= '{species}/{library_type}/deseq/dge.tsv.gz',
        hl=os.path.join(workflow.basedir, 'data/highlight_genes.tsv'),
    output:
        volcano= '{species}/{library_type}/deseq/volcano.pdf',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)
library(ggrepel)

dge <- fread('{input.dge}')
hl <- fread('{input.hl}')
dge[, highlight := gene_id %in% hl$gene_id]
dge[, contrast := gsub('_', ' ', contrast)]

gg <- ggplot(data=dge, aes(x=log2FoldChange, y=-log10(pvalue))) +
    geom_vline(xintercept=0, colour='dodgerblue', linetype='dashed', linewidth=0.25) +
    geom_point(data=dge[padj >= 0.01, ], pch='.', colour='grey80') +
    geom_point(data=dge[padj < 0.01, ], pch='.', colour='black') +
    geom_point(data=dge[highlight == TRUE & padj < 0.01, ], pch=19, size=0.75, colour='dodgerblue') +
    geom_text_repel(data=dge[highlight == TRUE & padj < 0.01,], aes(label=gene_name), colour='grey20', size=3) +
    scale_x_continuous(breaks=scales::pretty_breaks(10)) +
    xlab('log2 fold-change') +
    facet_wrap(~contrast) +
    theme_light() +
    theme(strip.text=element_text(colour='black'))
ggsave('{output.volcano}', width=12, height=10, units='cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


