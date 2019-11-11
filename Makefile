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

references/GRCh38.6:
	echo "Downloading snpEff database..."
	$(DOCKER_RUN) \
		quay.io/biocontainers/snpeff@sha256:5c61b86bf531d3bf20c0fe50e8197a35b977c281ef74369c67e842eb4d092941 \
		java -jar /usr/local/share/snpeff-4.3.1t-1/snpEff.jar download -dataDir /references GRCh38.86

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
		jmonlong/sveval@sha256:09d1ac8c942eca62a0e68385ac3624425d0a71004c73bf584aa9af9048847303 \
		R -e "sveval::freqAnnotate('/data/$(PREREQ)', '/references/gnomad_v2_sv.sites.pass.lifted.vcf.gz', out.vcf='/data/$(TARGET)', min.cov=.1)"

%.ann.vcf: %.vcf 
	echo "Annotating variants..."
	$(DOCKER_RUN) \
		quay.io/biocontainers/snpeff@sha256:5c61b86bf531d3bf20c0fe50e8197a35b977c281ef74369c67e842eb4d092941 \
		java -Xmx10000m -jar /usr/local/share/snpeff-4.3.1t-1/snpEff.jar \
		-dataDir /references \
		-t -quiet \
		-noNextProt -noMotif -noStats -classic \
		-no PROTEIN_PROTEIN_INTERACTION_LOCUS -no PROTEIN_STRUCTURAL_INTERACTION_LOCUS \
		GRCh38.86 /data/$(PREREQ) > $(@)
	chown `id -u`:`stat -c "%g" samples/` $(@)

#
# Reports
#

%.sv-report.pdf: %.sniffles.ann.freqGnomADcov10.vcf %.svim.ann.freqGnomADcov10.vcf references/gene_position_info.txt references/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz references/simpleRepeat.txt.gz sv-report.Rmd
	echo "Producing SV report..."
	$(DOCKER_RUN) \
		-v `realpath .`:/app -w /app \
		jmonlong/sveval-rmarkdown@sha256:78ec614cbdb46f64d76794ed2eda2a4f61fc6216dbb6a9d8843a3d7693ec8c1a \
		Rscript -e 'rmarkdown::render("sv-report.Rmd", output_format="pdf_document")' /data/$$(echo $^ | cut -f1 -d' ' | xargs basename) /data/$$(echo $^ | cut -f2 -d' ' | xargs basename) /references/gene_position_info.txt /references/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz /references/simpleRepeat.txt.gz
	mv sv-report.pdf $@

%.sv-report.html: %.sniffles.ann.freqGnomADcov10.vcf %.svim.ann.freqGnomADcov10.vcf references/gene_position_info.txt references/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz references/simpleRepeat.txt.gz
	echo "Producing SV report..."
	$(DOCKER_RUN) \
		-v `realpath .`:/app -w /app \
		jmonlong/sveval-rmarkdown@sha256:d2f504ee111aeaeecdb061d0adf86aca766a8ab116c541ce208b79bc9c448cbc \
		Rscript -e 'rmarkdown::render("sv-report.Rmd", output_format="html_document")' /data/$$(echo $^ | cut -f1 -d' ' | xargs basename) /data/$$(echo $^ | cut -f2 -d' ' | xargs basename) /references/gene_position_info.txt /references/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz /references/simpleRepeat.txt.gz
	mv sv-report.html $@
