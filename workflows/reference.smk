rule download_fasta_reference:
    output:
        genome= '{species}/ref/genome.fasta',
        fai= '{species}/ref/genome.fasta.fai',
    params:
        url= lambda wc: species[species.species== wc.species].fasta.iloc[0],
        gunzip= lambda wc: '| gunzip' if species[species.species== wc.species].fasta.iloc[0].endswith('.gz') else '',
    shell:
        r"""
        curl -s -L {params.url} {params.gunzip} > {output.genome}
        samtools faidx {output.genome}
        """


rule download_annotation:
    output:
        gxf= '{species}/ref/annotation.gtf',
    params:
        url= lambda wc: species[species.species== wc.species].gff.iloc[0],
        gunzip= lambda wc: '| gunzip' if species[species.species== wc.species].fasta.iloc[0].endswith('.gz') else '',
    shell:
        r"""
        curl -s -L {params.url} {params.gunzip} > {output.gxf}
        """


rule gene_description:
    input:
        gff= '{species}/ref/annotation.gff',
    output:
        descr= '{species}/ref/gene_description.tsv',
    shell:
        r"""
        {workflow.basedir}/scripts/getGffAttributes.py --gff {input.gff} -t gene -a ID Name description -r ID:gene_id > {output.descr}
        """
