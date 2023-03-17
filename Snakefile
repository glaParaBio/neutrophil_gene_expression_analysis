import pandas

def get_fastq(wc):
    fq = [ssfq[ssfq.fastq_base == wc.fastq_base].fastq_r1.iloc[0]]
    if not ssfq[ssfq.fastq_base == wc.fastq_base].fastq_r2.isna().values.any():
        fq2 = ssfq[ssfq.fastq_base == wc.fastq_base].fastq_r2.iloc[0]
        fq.append(fq2)
    return fq

species = pandas.read_csv(config['species'], sep= '\t', comment= '#')

ssfq = pandas.read_csv(config['sample_sheet'], sep= '\t', comment= '#', dtype= {'library_id': str, 'run_id': str})
ss = ssfq.drop(['run_id', 'fastq_r1', 'fastq_r2'], axis= 1).drop_duplicates()
assert len(ss.library_id) == len(set(ss.library_id))

for idx,row in ssfq.iterrows():
    for r12 in ['fastq_r1', 'fastq_r2']:
        if not row[r12].startswith('sra/'):
            ssfq.loc[idx, r12] = os.path.join(config['fastqdir'], row[r12])

ssfq['fastq_base'] = ssfq.library_id + '.' + ssfq.run_id 
assert len(ssfq.fastq_base) == len(set(ssfq.fastq_base))

wildcard_constraints:
    library_id= '|'.join([re.escape(x) for x in ss.library_id]),
    fastq_base= '|'.join([re.escape(x) for x in ssfq.fastq_base]),
    species=  '|'.join([re.escape(x) for x in ss.species.unique()]),


rule all:
    input:
        expand('{species}/{library_type}/bigwig/{library_id}.bw', zip, 
                library_type= ss[ss.library_type == 'RNAseq'].library_type, 
                species= ss[ss.library_type == 'RNAseq'].species, 
                library_id= ss[ss.library_type == 'RNAseq'].library_id),
        expand('{species}/{library_type}/featureCounts/{library_id}.mapq10.tsv', zip, 
                library_type= ss[ss.library_type == 'RNAseq'].library_type, 
                species= ss[ss.library_type == 'RNAseq'].species, 
                library_id= ss[ss.library_type == 'RNAseq'].library_id),
        'Homo_sapiens.GRCh38.105/RNAseq/clusterprofiler/gsea.tsv.gz',
        os.path.join(workflow.basedir, 'results/pca.pdf'),
        'Homo_sapiens.GRCh38.105/RNAseq/deseq/dge.tsv.gz',
        'Homo_sapiens.GRCh38.105/RNAseq/deseq/volcano.pdf',
        'Homo_sapiens.GRCh38.105/RNAseq/deseq/gex_vsd.tsv.gz',
        'Homo_sapiens.GRCh38.105/RNAseq/misc/hif1a_msigdb_targets_bees.pdf',

include: 'workflows/reference.smk'
include: 'workflows/qc.smk'
include: 'workflows/alignment.smk'
include: 'workflows/dge.smk'


rule prepare_fastq:
    input:
        r1= lambda wc: ssfq[ssfq.fastq_base == wc.fastq_base].fastq_r1.iloc[0],
        r2= lambda wc: ssfq[ssfq.fastq_base == wc.fastq_base].fastq_r2.iloc[0],
    output:
        r1= temp('fastq/{fastq_base}.R1.fastq.gz'),
        r2= temp('fastq/{fastq_base}.R2.fastq.gz'),
    params:
        pattern= lambda wc: ssfq[ssfq.fastq_base == wc.fastq_base].clontech_barcode.iloc[0],
    shell:
        # NB: The barcode is read 2 so we use R2 as primary input to umi_tools
        # and R1 as second
        r"""
        umi_tools extract --bc-pattern='{params.pattern}' -I {input.r2} -S {output.r2} \
            --read2-in {input.r1} --read2-out {output.r1}
        """
