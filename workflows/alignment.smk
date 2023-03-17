def featureCountsStrand(wc):
    """Numeric code used by featureCounts to identify strandness of protocol
    for RNAseq
    """
    strand= ss[ss.library_id == wc.library_id].strand
    assert len(strand) == 1
    strand= strand.iloc[0]
    code= None
    if strand == 'FR':
        code= 1
    elif strand == 'RF':
        code= 2
    elif strand == 'unstranded':
        code= 0
    else:
        sys.stderr.write('\n\nInvalid strand identifier found in sample sheet: %s\n' % strand)
        sys.exit(1)
    return code


rule hisat_index:
    input:
        species= '{species}/ref/genome.fasta',
    output:
        idx= '{species}/ref/genome.8.ht2',
    run:
        idx = re.sub('\.8\.ht2$', '', output.idx)

        shell(f"""
        hisat2-build -p 4 --seed 1234 -f {input.species} {idx}
        """)


rule hisat2:
    input:
        fq= ['fastq/{fastq_base}.R1.fastq.gz', 'fastq/{fastq_base}.R2.fastq.gz'],
        idx= lambda wc: '{species}/ref/genome.8.ht2',
    output:
        bam= temp('{species}/RNAseq/hisat2/{fastq_base}.bam'),
        hlog= '{species}/RNAseq/hisat2/{fastq_base}.log',
    params:
        max_intron_len= lambda wc: species[species.species == wc.species].max_intron_len.iloc[0] ,
    run:
        if len(input.fq) == 1:
            raw_fq = f'-U {input.fq[0]}'
        elif len(input.fq) == 2:
            raw_fq = f'-1 {input.fq[0]} -2 {input.fq[1]}'
        else:
            raise Exception('Unexpected number of fastq files')
        
        idx = re.sub('\.8\.ht2$', '', input.idx)

        shell(f"""
        hisat2 --summary-file {output.hlog} --new-summary --fr --rna-strandness RF \
           --max-intronlen {params.max_intron_len} --threads 4 -x {idx} {raw_fq} \
        | samtools view -u -@ 4 \
        | samtools sort -@ 8 > {output.bam}
        """)


rule merge_hisat2:
    input:
        bam= lambda wc: expand('{{species}}/RNAseq/hisat2/{fastq_base}.bam', fastq_base= ssfq[ssfq.library_id == wc.library_id].fastq_base),
    output:
        bam= temp('{species}/RNAseq/hisat2/{library_id}.dup.bam'),
    run:
        if len(input.bam) == 1:
            shell("mv {input.bam} {output.bam}")
        else:
            shell("samtools merge {output.bam} {input.bam}")


rule index_hisat_bam:
    input:
        bam= '{species}/RNAseq/hisat2/{library_id}.dup.bam',
    output:
        temp('{species}/RNAseq/hisat2/{library_id}.dup.bam.bai'),
    shell:
        'samtools index -@ 4 {input.bam}'


rule dedup:
    input:
        bam= '{species}/RNAseq/hisat2/{library_id}.dup.bam',
        bai= '{species}/RNAseq/hisat2/{library_id}.dup.bam.bai',
    output:
        bam= '{species}/RNAseq/hisat2/{library_id}.bam',
        xlog= '{species}/RNAseq/hisat2/{library_id}.dedup.log',
    shell:
        r"""
        umi_tools dedup --temp-dir . --paired -I {input.bam} -S {output.bam} -L {output.xlog} 
        """

rule dedup_index:
    input:
        bam= '{species}/RNAseq/hisat2/{library_id}.bam',
    output:
        bai= '{species}/RNAseq/hisat2/{library_id}.bam.bai',
    shell:
        r"""
        samtools index -@ 4 {input.bam}
        """


rule bigwig:
    input:
        bam= '{species}/RNAseq/hisat2/{library_id}.bam',
        bai= '{species}/RNAseq/hisat2/{library_id}.bam.bai',
    output:
        bw= '{species}/RNAseq/bigwig/{library_id}.bw',
    shell:
        r"""
        bamCoverage -b {input.bam} -o {output} \
            --binSize 50 \
            --minMappingQuality 5 \
            --normalizeUsing BPM \
            --numberOfProcessors 4
        """


rule count_reads_in_genes:
    input:
        bam= '{species}/RNAseq/hisat2/{library_id}.bam',       
        gff= '{species}/ref/annotation.gtf',
    output:
        counts= '{species}/RNAseq/featureCounts/{library_id}.mapq{mapq}.tsv',
    params:
        strand= featureCountsStrand,
        gene_id_key= 'gene_id',
        mapq= lambda wc: '-Q %s' % wc.mapq if int(wc.mapq) > 0 else '-Q 0 --primary -M'
    shell:
        r"""
        featureCounts -p -T 8 {params.mapq} --extraAttributes gene_name,gene_biotype \
            -s {params.strand} -t exon -g {params.gene_id_key} -a {input.gff} -o {output.counts} {input.bam}
        """


rule collate_counts:
    input:
        counts= lambda wc: [f'{wc.species}/RNAseq/featureCounts/{x}.tsv' for x in 
                sorted(ss[(ss.species == wc.species) & (ss.library_type == 'RNAseq')].library_id)],
    output:
        counts= '{species}/RNAseq/featureCounts/counts.tsv.gz',
    run:
        dts = []
        for x in input.counts:
            dt = pandas.read_csv(x, sep= '\t', comment= '#')
            fid = dt.columns[len(dt.columns)-1]
            dt.rename(columns= {fid: 'count'}, inplace= True)
            dt.rename(columns= {'Geneid': 'gene_id'}, inplace= True)
            library_id = re.sub('\.tsv$', '', os.path.basename(x))
            dt['library_id'] = library_id
            assert library_id in os.path.basename(fid)
            
            dts.append(dt[['library_id', 'gene_id', 'Length', 'count']])

        counts = pandas.concat(dts)
        counts.to_csv(output.counts, sep= '\t', index= False, na_rep= 'NA')

rule quality:
    input:
        ss= config['sample_sheet'],
        cnt= expand('{species}/{library_type}/featureCounts/{library_id}.mapq0.tsv', zip, 
                library_type= ss[ss.library_type == 'RNAseq'].library_type, 
                species= ss[ss.library_type == 'RNAseq'].species, 
                library_id= ss[ss.library_type == 'RNAseq'].library_id),
    output:
        pca= os.path.join(workflow.basedir, 'results/pca.pdf'),
        vsd= 'Homo_sapiens.GRCh38.105/RNAseq/deseq/gex_vsd.tsv.gz',
        biotypes= os.path.join(workflow.basedir, 'results/biotypes.tsv'),
        corr= os.path.join(workflow.basedir, 'results/correlations.pdf'),
        gex= os.path.join(workflow.basedir, 'results/histogram_gex.pdf'),
        count_summary= os.path.join(workflow.basedir, 'results/count_summary.tsv'),
        cnt= os.path.join(workflow.basedir, 'results/count_all_primary.tsv.gz'),
    script:
        '../scripts/quality.R'
