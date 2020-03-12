#
#	References:
#	https://monades.roperzh.com/rediscovering-make-automatic-variables/
#
# Keep all intermediate artifacts, should add a few to prune later...
.SECONDARY:

# Path to the Makefile
APP_PATH ?= $(dir $(realpath $(firstword $(MAKEFILE_LIST))))

# filename of first prerequisite ex. sample.fq
PREREQ = $(<F)

# filename of the target generated from the prerequisite ex. sample.bam
TARGET = $(@F)

CPU ?= 32

# Run as the calling user and samples directory group and
# map the sample into /data and references into /references
DOCKER_RUN = docker run -it --rm --cpus="$(CPU)" \
		--user `id -u`:`stat -c "%g" samples/` \
		-v `realpath references`:/references \
		-v `realpath $(@D)`:/data

#
# General recipes
#

%.md5: %
	echo "Calculating md5 of $(<) for provenance..."
	md5sum $(<) > $(<).md5

#
# Download references
#

references/hg38.fa:
	echo "Downloading hg38 reference and indexing..."
	mkdir -p references
	wget -N -P references https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	gunzip references/hg38.fa.gz
	md5sum -c $(APP_PATH)hg38.fa.md5

references/hg38.fa.fai: references/hg38.fa
	$(DOCKER_RUN) \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		faidx /references/hg38.fa

references/gene_position_info.txt:
	echo "Downloading gene list..."
	mkdir -p references
	wget -N -P references https://github.com/ucsc-upd/operations/files/3628601/gene_position_info.txt

references/gene_position_info.tsv:
	echo "/!\  Download it manually, I haven't put it anywhere yet /!\ "

references/gnomad_v2_sv.sites.pass.lifted.vcf.gz:
	echo "Downloading SV catalog from gnomAD-SV..."
	mkdir -p references
	wget -N -P references https://storage.googleapis.com/jmonlong-vg-wdl-dev-test/gnomad_v2_sv.sites.pass.lifted.vcf.gz

references/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz:
	echo "Downloading LoF intolerance score from gnomAD..."
	mkdir -p references
	wget -N -P references https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

references/simpleRepeat.txt.gz:
	echo "Downloading repeat annotation..."
	mkdir -p references
	wget -N -P references https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz

references/hsvlr.vcf.gz:
	echo "Downloading public SV catalog from long-read studies..."
	mkdir -p references
	wget -N -P references https://storage.cloud.google.com/jmonlong-vg-sv/hsvlrnodupins/hsvlrnodupins.vcf.gz
	mv references/hsvlrnodupins.vcf.gz references/hsvlr.vcf.gz

references/GRCh38_hg38_variants_2016-08-31.txt:
	echo "Downloading DGV SV catalog..."
	mkdir -p references
	wget -N -P references http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2016-08-31.txt

references/iscaPathogenic.txt.gz:
	echo "Downloading ClinGen Pathogenic CNV catalog from UCSC Genome Browser..."
	mkdir -p references
	wget -N -P references https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/iscaPathogenic.txt.gz

references/ENCFF010WHH.bed.gz:
	echo "Downloading CTCF peaks in kidney from ENCODE..."
	mkdir -p references
	wget -N -P references https://www.encodeproject.org/files/ENCFF010WHH/@@download/ENCFF010WHH.bed.gz

references/ENCFF166QIT.bed.gz:
	echo "Downloading regulatory region in kidney from ENCODE (hg19)..."
	mkdir -p references
	wget -N -P references https://www.encodeproject.org/files/ENCFF166QIT/@@download/ENCFF166QIT.bed.gz

references/ENCFF166QIT.lifted.bed.gz: references/ENCFF166QIT.bed.gz
	echo "Lift over regulatory regions to GCRh38..."
	mkdir -p references
	gunzip -f references/ENCFF166QIT.bed.gz
	$(DOCKER_RUN) \
		jmonlong/liftover@sha256:c1e513c7bede70edf4bb758edb061a457f4e76f2c0149993eca36baa486d235b \
		liftOver \
		/references/ENCFF166QIT.bed /home/hg19ToHg38.over.chain.gz \
		/references/ENCFF166QIT.lifted.bed /references/ENCFF166QIT.notlifted.bed
	gzip -f references/ENCFF166QIT.bed references/ENCFF166QIT.lifted.bed references/ENCFF166QIT.notlifted.bed

references/GRCh38.86:
	echo "Downloading snpEff database..."
	$(DOCKER_RUN) \
		quay.io/biocontainers/snpeff@sha256:5c61b86bf531d3bf20c0fe50e8197a35b977c281ef74369c67e842eb4d092941 \
		java -jar /usr/local/share/snpeff-4.3.1t-1/snpEff.jar download -dataDir /references GRCh38.86
	zcat references/GRCh38.86/regulation_DNase1-NHEK_enriched_sites.bin | sed "s/ /_/g" | gzip > references/GRCh38.86/regulation_DNase1-NHEK_enriched_sites_nospace.bin
	zcat references/GRCh38.86/regulation_DNase1-HUVEC_enriched_sites.bin | sed "s/ /_/g" | gzip > references/GRCh38.86/regulation_DNase1-HUVEC_enriched_sites_nospace.bin


references/GRCh38.86/phastCons: references/GRCh38.86 references/hg38.fa.fai
	echo "Downloading snpEff conservation database..."
	mkdir -p references/GRCh38.86/phastCons
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr1.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr2.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr3.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr4.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr5.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr6.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr7.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr8.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr9.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr10.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr11.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr12.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr13.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr14.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr15.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr16.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr17.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr18.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr19.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr20.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr21.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr22.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chrM.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chrX.phastCons100way.wigFix.gz
	wget -N -P references/GRCh38.86/phastCons http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chrY.phastCons100way.wigFix.gz
	cp references/hg38.fa.fai references/GRCh38.86/phastCons/genome.fai


#
# Download NA12878 chr11 from https://github.com/nanopore-wgs-consortium
# and convert to fq as a test sample
#
# For a whole genome version download:
# http://s3.amazonaws.com/nanopore-human-wgs/rel6/rel_6.fastq.gz
#

samples/na12878-chr11/na12878-chr11.fq.gz:
	echo "Downloading na12878 chr11 bam..."
	mkdir -p samples/na12878-chr11 references
	wget -O samples/na12878-chr11/na12878-chr11.original.bam \
		http://s3.amazonaws.com/nanopore-human-wgs/chr11.sorted.bam
	echo "Converting bam to fq..."
	$(DOCKER_RUN) \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		fastq na12878-chr11.original.bam -0 na12878-chr11.fq.gz

#
# Map sample to hg38 creating a sorted and indexed bam
#

%.sam: %.fq.gz %.fq.gz.md5 references/hg38.fa references/hg38.fa.fai
	echo "Mapping reads to reference genome..."
	$(DOCKER_RUN) \
		tpesout/minimap2@sha256:5df3218ae2afebfc06189daf7433f1ade15d7cf77d23e7351f210a098eb57858 \
		-ax map-ont --MD -t $(CPU) /references/hg38.fa $(PREREQ)
	mv $(@D)/minimap2.sam $(@)

%.bam: %.sam
	$(DOCKER_RUN) \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		view -S -b $(PREREQ) -o $(TARGET) -@ $(CPU)

%.sorted.bam: %.bam
	$(DOCKER_RUN) \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		sort $(PREREQ) -o $(TARGET) -@ $(CPU)

%.sorted.bam.bai: %.sorted.bam
	$(DOCKER_RUN) \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		index $(PREREQ)

#
# Call variants against hg38
#

%.sniffles.vcf: %.sorted.bam references/hg38.fa
	echo "Calling variants with sniffles..."
	$(DOCKER_RUN) \
		quay.io/biocontainers/sniffles@sha256:98a5b91db2762ed3b8aca3bffd3dca7d6b358d2a4a4a6ce7579213d9585ba08a \
		sniffles -m /data/$(PREREQ) -v /data/$(TARGET) --genotype

%.svim.vcf: %.sorted.bam %.sorted.bam.bai references/hg38.fa
	echo "Calling variants with svim..."
	$(DOCKER_RUN) \
		quay.io/biocontainers/svim@sha256:4239718261caf12f6c27d36d5657c13a2ca3b042c833058d345b04531b576442 \
		svim alignment /data/svim /data/$(PREREQ) /references/hg38.fa --sample $(TARGET)
	mv $(@D)/svim/final_results.vcf $(@)

%.freqGnomADcov10.vcf: %.vcf references/gnomad_v2_sv.sites.pass.lifted.vcf.gz
	echo "Annotating SV frequency using the gnomAD-SV catalog..."
	$(DOCKER_RUN) \
		jmonlong/sveval@sha256:fe35f70ae2765bddaf0a990d9f04f94b029c7006ebabab653d046da4566f3fbb \
		R -e "sveval::freqAnnotate('/data/$(PREREQ)', '/references/gnomad_v2_sv.sites.pass.lifted.vcf.gz', out.vcf='/data/$(TARGET)', min.cov=.1)"

%.ann.vcf: %.vcf references/GRCh38.86
	echo "Annotating variants..."
	$(DOCKER_RUN) \
		quay.io/biocontainers/snpeff@sha256:5c61b86bf531d3bf20c0fe50e8197a35b977c281ef74369c67e842eb4d092941 \
		java -Xmx10000m -jar /usr/local/share/snpeff-4.3.1t-1/snpEff.jar \
		-dataDir /references \
		-t -quiet \
		-noNextProt -noMotif -noStats -classic \
		-reg DNase1-HUVEC_enriched_sites_nospace -reg DNase1-NHEK_enriched_sites_nospace \
		-no PROTEIN_PROTEIN_INTERACTION_LOCUS -no PROTEIN_STRUCTURAL_INTERACTION_LOCUS \
		GRCh38.86 /data/$(PREREQ) > $(@)
	chown `id -u`:`stat -c "%g" samples/` $(@)

%.sift.vcf: %.vcf references/GRCh38.86/phastCons
	echo "Annotating variants..."
	$(DOCKER_RUN) \
		quay.io/biocontainers/snpsift@sha256:2b2a0fa662bde7bd8643191c02eb35a3abd7bb426928b722bcfc4edd4e66d87d \
		java -Xmx10000m -jar /usr/local/share/snpsift-4.2-4/SnpSift.jar \
		phastCons /references/GRCh38.86/phastCons /data/$(PREREQ) > $(@)
	chown `id -u`:`stat -c "%g" samples/` $(@)

#
# Reports
#

%.sv-report.html: %.sniffles.ann.freqGnomADcov10.vcf %.svim.ann.freqGnomADcov10.vcf references/gene_position_info.tsv references/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz references/simpleRepeat.txt.gz references/hsvlr.vcf.gz references/GRCh38_hg38_variants_2016-08-31.txt references/iscaPathogenic.txt.gz references/ENCFF010WHH.bed.gz references/ references/gnomad_v2_sv.sites.pass.lifted.vcf.gz references/ENCFF166QIT.lifted.bed.gz
	echo "Producing SV report..."
	$(DOCKER_RUN) \
		-v `realpath .`:/app -w /app \
		jmonlong/sveval-rmarkdown@sha256:99e92947e226ae8496e9b061701ff5d89a445b61c919bd26744756d6b97d6a69 \
		Rscript -e 'rmarkdown::render("sv-report.Rmd", output_format="html_document")' /data/$$(echo $^ | cut -f1 -d' ' | xargs basename) /data/$$(echo $^ | cut -f2 -d' ' | xargs basename) /references/gene_position_info.tsv /references/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz /references/simpleRepeat.txt.gz /references/hsvlr.vcf.gz /references/GRCh38_hg38_variants_2016-08-31.txt /references/iscaPathogenic.txt.gz /references/ENCFF010WHH.bed.gz /references/gnomad_v2_sv.sites.pass.lifted.vcf.gz /references/ENCFF166QIT.lifted.bed.gz
	mv sv-report.html $@
	mv sv-te-like-insertions.tsv $@.sv-te-like-insertions.tsv
