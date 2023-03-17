base_fq= {}
for i, row in ssfq.iterrows():
    base_fq[os.path.basename(re.sub('\.fastq\.gz$|\.fq\.gz$', '', row['fastq_r1']))]= row['fastq_r1']
    if not pandas.isna(row['fastq_r2']):
        base_fq[os.path.basename(re.sub('\.fastq\.gz$|\.fq\.gz$', '', row['fastq_r2']))]= row['fastq_r2']


rule fastqc:
    priority: -10
    input:
        lambda wc: base_fq[wc.base_fq],
    output:
        'fastqc/{base_fq}_fastqc.zip',
    shell:
        r"""
        fastqc -o fastqc {input}
        """


rule multiqc_fastqc:
    priority: -10
    input:
        expand('fastqc/{base_fq}_fastqc.zip', base_fq= base_fq.keys()),
    output:
        'multiqc/fastqc_report.html',
    shell:
        r"""
        multiqc --force --outdir `dirname {output}` \
            --filename `basename {output}` {input}
        """
