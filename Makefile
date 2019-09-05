# Sample ID - will be used to prefix all files generated
ID ?= test

# Max number of CPUs/Cores to use
CPU ?= 32

APP_PATH ?= $(dir $(realpath $(firstword $(MAKEFILE_LIST))))

#
# Download references, test samples, and generate MD5's...
#

data/references/hg38.fa:
	echo "Downloading hg38 reference and indexing..."
	mkdir -p data/references
	wget -N -P data/references https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	gunzip data/references/hg38.fa.gz
	md5sum -c $(APP_PATH)hg38.fa.md5

data/references/hg38.fa.fai: data/references/hg38.fa
	docker run -it --rm --cpus="$(CPU)" \
		-v `realpath data/references`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		faidx /data/hg38.fa

data/$(ID)/$(ID).fq data/$(ID)/$(ID).fa:
	echo "Downloading GM24385.chr20.fq reference quality nanopore fastq as a test sample..."
	mkdir -p data/$(ID)
	wget -O data/$(ID)/$(ID).fq https://lc2019.s3-us-west-2.amazonaws.com/sample_data/GM24385/GM24385.chr20.fq
	md5sum -c $(APP_PATH)test.fq.md5
	sed -n '1~4s/^@/>/p;2~4p' data/$(ID)/$(ID).fq > data/$(ID)/$(ID).fa

data/$(ID)/$(ID).fq.md5: data/$(ID)/$(ID).fq
	echo "Calculating MD5 of FASTQ for provenance..."
	md5sum data/$(ID)/$(ID).fq > data/$(ID)/$(ID).fq.md5

#
# Map sample to hg38
#

data/$(ID)/$(ID).minimap2_hg38.sam: \
	data/references/hg38.fa data/references/hg38.fa.fai data/$(ID)/$(ID).fq data/$(ID)/$(ID).fq.md5
	echo "Mapping reads to reference genome..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		tpesout/minimap2@sha256:5df3218ae2afebfc06189daf7433f1ade15d7cf77d23e7351f210a098eb57858 \
		-ax map-ont -t $(CPU) /references/hg38.fa /data/$(ID).fq
	mv data/$(ID)/minimap2.sam data/$(ID)/$(ID).minimap2_hg38.sam
	mv data/$(ID)/minimap2.log data/$(ID)/$(ID).minimap2_hg38.log

data/$(ID)/$(ID).minimap2_hg38_sorted.bam: data/$(ID)/$(ID).minimap2_hg38.sam
	echo "Converting sam to sorted bam with index..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		view -S -b /data/$(ID).minimap2_hg38.sam -o /data/$(ID).minimap2_hg38.bam
	echo "Sorting bam..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		sort /data/$(ID).minimap2_hg38.bam -o /data/$(ID).minimap2_hg38_sorted.bam
	echo "Indexing bam..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		index /data/$(ID).minimap2_hg38_sorted.bam

#
# Call variants against hg38
#

data/references/trainedModels:
	cd data/references && curl http://www.bio8.cs.hku.hk/trainedModels.tbz | tar -jxf -

data/$(ID)/$(ID).hg38_lite.vcf: \
	data/$(ID)/$(ID).minimap2_hg38_sorted.bam data/references/hg38.fa data/references/trainedModels
	echo "Calling variants against hg38..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		--user=`id -u`:`id -g` \
		-e CUDA_VISIBLE_DEVICES="" \
		lifebitai/clairvoyante:latest \
		pypy /opt/conda/bin/clairvoyante/callVarBam.py \
			--chkpnt_fn /references/trainedModels/fullv3-illumina-novoalign-hg001+hg002-hg38/learningRate1e-3.epoch500 \
			--ref_fn /references/hg38.fa \
			--bam_fn /data/$(ID).minimap2_hg38_sorted.bam \
			--call_fn /data/$(ID).hg38_lite.vcf \
			--ctgName chr20 \
			--ctgStart 1000000 \
			--ctgEnd 1010000

benchmark:
	# Compare the region we tried to call variants on to the ground truth from GIAB
	# First download the ground truth
	wget -N -P data/references https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
	wget -N -P data/references https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi
	wget -N -P data/references ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
	# Compress and index our variant set - should be up in variant calling?
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		--user=`id -u`:`id -g` \
		quay.io/biocontainers/bcftools@sha256:d9360c33b0771a49f49db174d03e0e49c8720eef34330be4c0371ca54fbb8f7b \
		bcftools view -Oz -o /data/$(ID).hg38_lite.vcf.gz /data/$(ID).hg38_lite.vcf
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		--user=`id -u`:`id -g` \
		quay.io/biocontainers/bcftools@sha256:d9360c33b0771a49f49db174d03e0e49c8720eef34330be4c0371ca54fbb8f7b \
		bcftools index /data/$(ID).hg38_lite.vcf.gz
	# Find the variants found in both
	mkdir -p data/$(ID)/diff
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		--user=`id -u`:`id -g` \
		quay.io/biocontainers/bcftools@sha256:d9360c33b0771a49f49db174d03e0e49c8720eef34330be4c0371ca54fbb8f7b \
		bcftools isec --regions chr20:1000000-1010000 \
		-p /data/diff -n=2 \
		/data/$(ID).hg38_lite.vcf.gz \
		/references/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
