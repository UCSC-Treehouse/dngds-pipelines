import os
from Bio import SeqIO

REF_DIR = config['ref_root']

##
## handy rules to run specific analysis
##

## CNV calling from Illumina reads using Control-FREEC
rule ill_cnvs:
    input: '{root}/freec/rd_baf_{w}bp/{samp}_{w}bp_CNVs'.format(root=config['root_ill'], samp=config['sample_ill'], w=config['bin_size'])

##
## SV report
##

if config['sample_ill'] != '':
    ## if we have illumina data, use the SV/CNV calls from it
    rule sv_report_cfg:
        input:
            sniffles_vcf='{root}/{sample}.sniffles.vcf',
            svim_vcf='{root}/{sample}.svim.vcf',
            cnv_ill='{root}/{sample}.freec_{w}bp.vcf'.format(w=config['bin_size'], sample=config['sample_ill'], root=config['root_ill']),
            smoove_ill='{root}/{sample}.smoove.vcf'.format(sample=config['sample_ill'], root=config['root_ill']),
            idxcov='{root}/indexcov_{sample}/indexcov_{sample}-indexcov.bed.gz',
            gene_pos= REF_DIR + '/gene_position_info.tsv',
            pli_gene= REF_DIR + '/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz',
            simprep= REF_DIR + '/simpleRepeat.txt.gz',
            hsvlr= REF_DIR + '/hsvlr.vcf.gz',
            dgv= REF_DIR + '/GRCh38_hg38_variants_2016-08-31.txt',
            clingen= REF_DIR + '/iscaPathogenic.txt.gz',
            ctcf= REF_DIR + '/ENCFF010WHH.bed.gz',
            gnomadsv= REF_DIR + '/gnomad_v2_sv.sites.pass.lifted.vcf.gz',
            cre= REF_DIR + '/ENCFF166QIT.lifted.bed.gz',
            cons= REF_DIR + '/phastConsElements100way.txt.gz',
            gencode= REF_DIR + '/gencode.v35.annotation.gtf.gz',
            cytoband= REF_DIR + '/cytoBandIdeo.txt.gz'
        output: '{root}/sv-report-{sample}.yaml'
        run:
            # prepare YAML config file
            cfg = ''
            for input_key in input.keys():
                cfg += '{}: "{}"\n'.format(input_key, input[input_key])
            with open(output[0], 'w') as cfg_outf:
                cfg_outf.write(cfg)
else:
    ## otherwise, just use LR calls
    rule sv_report_cfg:
        input:
            sniffles_vcf='{root}/{sample}.sniffles.vcf',
            svim_vcf='{root}/{sample}.svim.vcf',
            idxcov='{root}/indexcov_{sample}/indexcov_{sample}-indexcov.bed.gz',
            gene_pos= REF_DIR + '/gene_position_info.tsv',
            pli_gene= REF_DIR + '/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz',
            simprep= REF_DIR + '/simpleRepeat.txt.gz',
            hsvlr= REF_DIR + '/hsvlr.vcf.gz',
            dgv= REF_DIR + '/GRCh38_hg38_variants_2016-08-31.txt',
            clingen= REF_DIR + '/iscaPathogenic.txt.gz',
            ctcf= REF_DIR + '/ENCFF010WHH.bed.gz',
            gnomadsv= REF_DIR + '/gnomad_v2_sv.sites.pass.lifted.vcf.gz',
            cre= REF_DIR + '/ENCFF166QIT.lifted.bed.gz',
            cons= REF_DIR + '/phastConsElements100way.txt.gz',
            gencode= REF_DIR + '/gencode.v35.annotation.gtf.gz',
            cytoband= REF_DIR + '/cytoBandIdeo.txt.gz'
        output: '{root}/sv-report-{sample}.yaml'
        run:
            # prepare YAML config file
            cfg = ''
            for input_key in input.keys():
                cfg += '{}: "{}"\n'.format(input_key, input[input_key])
            with open(output[0], 'w') as cfg_outf:
                cfg_outf.write(cfg)   

rule sv_report:
    input: '{root}/sv-report-{sample}.yaml'
    output:
        html='{root}/sv-report-{sample}.html',
        tsv='{root}/sv-report-{sample}.tsv'
    singularity: "docker://jmonlong/sveval-rmarkdown@sha256:8cbd7de092884477f0620509b95420fedc226448e8e84dc62230fa728722b028"
    benchmark: '{root}/benchmarks/{sample}.svreport.tsv'
    log: '{root}/logs/{sample}.svreport.log'
    shell:
        """
        Rscript -e 'rmarkdown::render("sv-report.Rmd", output_format="html_document")' {input} {output.tsv} 2> {log}
        mv sv-report.html {output.html}
        """

rule unzip_vcf:
    input: '{vcf}.vcf.gz'
    output: '{vcf}.vcf'
    shell: "gunzip -c {input} > {output}"

##
## References
##

rule gencode_dwl:
    output: REF_DIR + '/gencode.v35.annotation.gtf.gz'
    shell:
        "wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz -O {output}"

rule cytoband_dwl:
    output: REF_DIR + '/cytoBandIdeo.txt.gz'
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz -O {output}"

rule hg38_dwl:
    output: REF_DIR + '/hg38.fa'
    params:
        temp_fa='hg38.fa.gz'
    shell:
        """
	wget -O {params.temp_fa} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	gunzip -c {params.temp_fa} > {output.fa}
        rm {params.temp_fa}
        """

rule index_fasta:
    input: '{file}.fa'
    output: '{file}.fa.fai'
    singularity: 'docker://quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f'
    shell: "samtools faidx {input}"

rule gnomadsv_dwl:
    output: REF_DIR + '/gnomad_v2_sv.sites.pass.lifted.vcf.gz'
    shell: "wget -O {output} https://storage.googleapis.com/jmonlong-vg-wdl-dev-test/gnomad_v2_sv.sites.pass.lifted.vcf.gz"

rule gnomad_pli_dwl:
    output: REF_DIR + '/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz'
    shell: "wget -O {output} https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

rule simprep_dwl:
    output: REF_DIR + '/simpleRepeat.txt.gz'
    shell: "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz"

rule hsvlr_dwl:
    output: REF_DIR + '/hsvlr.vcf.gz'
    shell: "wget -O {output} https://storage.cloud.google.com/jmonlong-vg-sv/hsvlrnodupins/hsvlrnodupins.vcf.gz"

rule dgv_dwl:
    output: REF_DIR + '/GRCh38_hg38_variants_2016-08-31.txt'
    shell: "wget -O {output} http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2016-08-31.txt"

rule clingen_dwl:
    output: REF_DIR + '/iscaPathogenic.txt.gz'
    shell: "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/iscaPathogenic.txt.gz"

rule ctcf_del:
    output: REF_DIR + '/ENCFF010WHH.bed.gz'
    shell: "wget -O {output} https://www.encodeproject.org/files/ENCFF010WHH/@@download/ENCFF010WHH.bed.gz"

rule cres_dwl:
    output: REF_DIR + '/ENCFF166QIT.bed.gz'
    shell: "wget -O {output} https://www.encodeproject.org/files/ENCFF166QIT/@@download/ENCFF166QIT.bed.gz"

rule chain_dwl:
    output: REF_DIR + '/hg19ToHg38.over.chain.gz'
    shell: "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

rule lift_over_bed:
    input:
           bed=REF_DIR + '/{region}.bed.gz',
           chain=REF_DIR + '/hg19ToHg38.over.chain.gz'
    output: REF_DIR + '/{region}.lifted.bed.gz'
    singularity: 'docker://jmonlong/liftover@sha256:c1e513c7bede70edf4bb758edb061a457f4e76f2c0149993eca36baa486d235b'
    params:
        temp_in='{region}.bed',
        temp_out='{region}.lifted.bed'
    shell:
        """
	gunzip -c {input.bed} > {params.temp_in}
	liftOver {params.temp_in}  {input.chain} {params.temp_out} notlifted.bed
	gzip -c {params.temp_out} > {output}
        rm -f {params.temp_in} {params.temp_out} notlifted.bed
        """

rule conselt_dwl:
    output: REF_DIR + '/phastConsElements100way.txt.gz'
    shell: "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/phastConsElements100way.txt.gz"

##
## Quick coverage information from goleft 
##

rule goleft_indexcov:
    input:
        bam='{root}/{sample}.sorted.bam',
        bai='{root}/{sample}.sorted.bam.bai'
    output: '{root}/indexcov_{sample}/indexcov_{sample}-indexcov.bed.gz'
    params:
        dir='{root}/indexcov_{sample}/'
    singularity: "docker://quay.io/biocontainers/goleft:0.2.0--0"
    benchmark: '{root}/benchmarks/{sample}.indexcov.tsv'
    log: '{root}/logs/{sample}.indexcov.log'
    shell: "goleft indexcov --directory {params.dir} --sex chrX,chrY {input.bam} 2> {log}"

##
## SV calling from nanopore data
##

rule call_sv_sniffles:
    input: '{root}/{sample}.sorted.bam'
    output: '{root}/{sample}.sniffles.vcf'
    singularity: 'docker://quay.io/biocontainers/sniffles@sha256:98a5b91db2762ed3b8aca3bffd3dca7d6b358d2a4a4a6ce7579213d9585ba08a'
    benchmark: '{root}/benchmarks/{sample}.sniffles.tsv'
    log: '{root}/logs/{sample}.sniffles.log'
    shell: "sniffles -m {input} -s 3 -v {output} --genotype 2> {log}"

rule call_sv_svim:
    input:
        bam='{root}/{sample}.sorted.bam',
        bai='{root}/{sample}.sorted.bam.bai',
        ref_fa='{}/{}.fa'.format(REF_DIR, config['ref']),
        ref_fai='{}/{}.fa.fai'.format(REF_DIR, config['ref'])
    output: '{root}/{sample}.svim.vcf'
    singularity: 'docker://quay.io/biocontainers/svim@sha256:7ae8dfc3fe9cce45aaa15ab56fe4ee93feee654f76219cec605709b70a1d44c2'
    benchmark: '{root}/benchmarks/{sample}.svim.tsv'
    log: '{root}/logs/{sample}.svim.log'
    params:
        res_dir='{root}/svim_results'
    shell:
        """
	svim alignment {params.res_dir} {input.bam} {input.ref_fa} --sample {wildcards.sample} 2> {log}
	mv {params.res_dir}/final_results.vcf {output}
        """

##
## CNV from Illumina reads using Control-FREEC
##
CHUNKS = range(1, config['nb_chunks']+1)
NBCHUNKS = config['nb_chunks']

rule cfg_freec_rd_baf:
    input:
        pileup = '{}/freec/{}.pileup.gz'.format(config['root_ill'],
                                                config['sample_ill']),
        bam=config['bam_ill'],
        chr_len = '{}/{}_chrs/chrs.len'.format(REF_DIR, config['ref']),
        snps = '{}/{}_snp{}.txt.gz'.format(REF_DIR,
                                           config['ref'],
                                           config['dbsnp_version'])
    output: '{}/freec/rd_baf_{{w}}bp/{}_{{w}}bp.cfg.txt'.format(config['root_ill'], config['sample_ill'])
    params:
        out_dir='{}/freec/rd_baf_{{w}}bp/'.format(config['root_ill']),
        chrs_dir='{}/{}_chrs/'.format(REF_DIR, config['ref'])
    run:
        # write Control-FREEC config file
        cfg = '[general]\nchrLenFile={chr_len}\nploidy=2\nbreakPointThreshold=.8\n'
        cfg += 'window={w}\nstep={w}\nchrFiles={chrs_dir}\nsex={sex}\noutputDir={outdir}\n'
        cfg += '[sample]\nmateFile={bam}\nminiPileup={pileup}\n'
        cfg += 'inputFormat=BAM\nmateOrientation=0\n[BAF]\nSNPfile={snps}\n'
        cfg = cfg.format(chr_len=input['chr_len'], w=wildcards.w,
                         chrs_dir=params.chrs_dir, sex=config['sex'],
                         outdir=params.out_dir,
                         bam=input['bam'], pileup=input['pileup'], snps=input['snps'])
        with open(output[0], 'w') as cfg_outf:
            cfg_outf.write(cfg)

rule freec_rd_baf:
    input:
        pileup = '{}/freec/{}.pileup.gz'.format(config['root_ill'],
                                                config['sample_ill']),
        bam = config['bam_ill'],
        chr_len = '{}/{}_chrs/chrs.len'.format(REF_DIR, config['ref']),
        snps = '{}/{}_snp{}.txt.gz'.format(REF_DIR,
                                           config['ref'],
                                           config['dbsnp_version']),
        cfg='{}/freec/rd_baf_{{w}}bp/{}_{{w}}bp.cfg.txt'.format(config['root_ill'], config['sample_ill'])
    output: '{}/freec/rd_baf_{{w}}bp/{}_{{w}}bp_CNVs'.format(config['root_ill'], config['sample_ill'])
    singularity: "docker://kfdrc/controlfreec:11.5"
    benchmark: '{}/benchmarks/{}.freec_{{w}}bp.tsv'.format(config['root_ill'], config['sample_ill'])
    log: '{}/logs/{}.freec_{{w}}bp.log'.format(config['root_ill'], config['sample_ill'])
    params:
        out_dir='{}/freec/rd_baf_{{w}}bp'.format(config['root_ill']),
        out_prefix=os.path.basename(config['bam_ill'])
    shell:
        """
        /FREEC-11.5/src/freec -conf {input.cfg} 2> {log}
        cp {params.out_dir}/{params.out_prefix}_CNVs {output}
        """

rule convert_freec_to_vcf:
    input:
        cnv='{root_ill}/freec/rd_baf_{w}bp/{sample}_{w}bp_CNVs',
        ref='{}/{}.fa'.format(REF_DIR, config['ref'])
    output: '{root_ill}/{sample}.freec_{w}bp.vcf'
    shell:
        "python3 convert_freec_to_vcf.py -v {input.cnv} -o {output} -r {input.ref}"

rule dwl_dbsnp:
    output: REF_DIR + '/{ref}_snp{v}.txt.gz'
    params: temp_txt=REF_DIR + '/{ref}_snp{v}.txt.temp.gz'
    shell:
        """
        wget -O {params.temp_txt} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.ref}/database/snp{wildcards.v}Common.txt.gz
        zcat {params.temp_txt} | awk 'BEGIN{{OFS="\\t"}}{{if(index($2, "_")==0 && $12=="single" && length($10)==3){{print $2,$3+1,$10,$7,$8,$5}}}}' | sort -k 1,1 -k 2,2n | gzip > {output}
        rm {params.temp_txt}
        """

rule snps_bed:
    input: REF_DIR + '/{ref}_snp{v}.txt.gz'
    output: REF_DIR + '/{ref}_snp{v}.bed'
    shell:
        """
        zcat {input} | awk 'BEGIN{{OFS="\t"}}{{print $1,$2-1,$2}}' > {output}
        """

rule chr_lens:
    input: '{root}/{ref}.fa'
    output: '{root}/{ref}_chrs/chrs.len'
    run:
        outf = open(output[0], 'w')
        for rec in SeqIO.parse(input[0], "fasta"):
            outf.write(rec.id + '\t' + len(rec.seq) + '\n')
            SeqIO.write([rec], os.path.dirname(output[0]) + '/' + rec.id + '.fa', 'fasta')
        outf.close()
        
rule merge_pileup:
    input: expand('{root}/freec/pileups/{label}.{chunk}.pileup', root=config['root_ill'], label=config['sample_ill'], chunk=CHUNKS)
    output: '{}/freec/{}.pileup.gz'.format(config['root_ill'], config['sample_ill'])
    shell:
        """
        cat {input} | sort -k 1,1 -k 2,2n | gzip > {output}
        """

rule pileup_chunk:
    input:
        bed='{}/{}_snp{}.bed'.format(REF_DIR, config['ref'], config['dbsnp_version']),
        bam=config['bam_ill'],
        ref='{}/{}.fa'.format(REF_DIR, config['ref'])
    output: '{}/freec/pileups/{}.{{chunk, \d+}}.pileup'.format(config['root_ill'], config['sample_ill'])
    container: "docker://quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f"
    params: temp_bed=config['root_ill'] + '/freec/pileups/temp_pileup_{chunk}.bed'
    benchmark: config['root_ill'] + '/freec/benchmark/pileup.' + config['sample_ill'] + '{chunk}.tsv'
    shell:
        """
        split -n l/{wildcards.chunk}/{NBCHUNKS} {input.bed} > {params.temp_bed}
        samtools mpileup -l {params.temp_bed} {input.bam} -f {input.ref} -Q 0 -q 1 > {output}
        rm {params.temp_bed}
        """
