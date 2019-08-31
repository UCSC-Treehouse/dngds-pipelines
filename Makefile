# Sample ID - will be used to prefix all files generated
ID ?= test

# Max number of CPUs/Cores to use
CPU ?= 32

# reserve some CPUs while running call_consensus
# CPU = 64, use 56, CPU = 32, use 28
HELEN_CALL_CONSENSUS_CPU ?= 28

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
	# Running but not working attempt to compare the lite vcf to a genome in a bottle variant set
	wget -N -P data/references https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
	wget -N -P data/references https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		--user=`id -u`:`id -g` \
		quay.io/biocontainers/tabix@sha256:b3c27ce674200727cd1b53a93ac4299f179e5411cd9c5a22db84b0c21b21baca \
		tabix -p vcf /references/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz \
		chr20:1000000-1010000 > data/references/HG002_GRCh38_GIAB_lite.vcf
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		--user=`id -u`:`id -g` \
		quay.io/biocontainers/snpsift@sha256:57ccff2c3f75f61990d7571c1d0e46255128bc118e39fd5a61c3fd84a440d1c9 \
		java -jar /usr/local/share/snpsift-4.2-3/SnpSift.jar concordance -v /references/HG002_GRCh38_GIAB_lite.vcf /data/$(ID).clairvoyante_hg38_lite.vcf

#
# De novo assemble and polish the sample
#

data/$(ID)/$(ID).shasta.fa: data/$(ID)/$(ID).fa data/$(ID)/$(ID).fq.md5
	echo "Assembling fasta..."
	rm -rf data/$(ID)/ShastaRun
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		tpesout/shasta@sha256:048f180184cfce647a491f26822f633be5de4d033f894ce7bc01e8225e846236 \
		--input /data/$(ID).fa
	mv data/$(ID)/ShastaRun/Assembly.fasta data/$(ID)/$(ID).shasta.fa

data/$(ID)/$(ID).minimap2_shasta.sam: data/$(ID)/$(ID).shasta.fa
	echo "Mapping fastq back to the assembled fasta..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		tpesout/minimap2@sha256:5df3218ae2afebfc06189daf7433f1ade15d7cf77d23e7351f210a098eb57858 \
		-ax map-ont -t $(CPU) /data/$(ID).shasta.fa /data/$(ID).fq
	mv data/$(ID)/minimap2.sam data/$(ID)/$(ID).minimap2_shasta.sam

data/$(ID)/$(ID).minimap2_shasta_sorted.bam: data/$(ID)/$(ID).minimap2_shasta.sam
	echo "Converting sam to sorted bam with index..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		view -S -b /data/$(ID).minimap2_shasta.sam -o /data/$(ID).minimap2_shasta.bam
	echo "Sorting bam..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		sort /data/$(ID).minimap2_shasta.bam -o /data/$(ID).minimap2_shasta_sorted.bam
	echo "Indexing bam..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		index /data/$(ID).minimap2_shasta_sorted.bam

data/$(ID)/$(ID).marginPolish.fa: data/$(ID)/$(ID).minimap2_shasta_sorted.bam data/$(ID)/$(ID).shasta.fa
	echo "Polishing fasta..."
	rm -rf data/$(ID)/marginPolish
	mkdir -p data/$(ID)/marginPolish
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		tpesout/margin_polish@sha256:de10c726bcc6af2f58cbb35af32ed0f0d95a3dc5f64f66dcc4eecbeb36f98b65 \
		/data/$(ID).minimap2_shasta_sorted.bam /data/$(ID).shasta.fa \
		/opt/MarginPolish/params/allParams.np.human.guppy-ff-235.json -t $(CPU) -o /data/marginPolish/output -f
	mv data/$(ID)/marginPolish/output.fa data/$(ID)/$(ID).marginPolish.fa

data/$(ID)/$(ID).assembly.fa: data/$(ID)/$(ID).marginPolish.fa
	echo "Downloading Helen model..."
	wget -N -P data/references https://storage.googleapis.com/kishwar-helen/helen_trained_models/v0.0.1/r941_flip235_v001.pkl
	echo "Calling consensus via Helen..."
	rm -rf helen_hdf5/
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		kishwars/helen@sha256:ac3504a6450c57b138a51652ebfff39cf1f5a0435ec995fcf137c7a1978cbedd \
		call_consensus.py \
		-i /data/marginPolish/ \
		-m /references/r941_flip235_v001.pkl \
		-o helen_hdf5/ \
		-p prediction \
		-w 0 \
		-t $(HELEN_CALL_CONSENSUS_CPU)
	echo "Stitching via Helen..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		kishwars/helen@sha256:ac3504a6450c57b138a51652ebfff39cf1f5a0435ec995fcf137c7a1978cbedd \
		stitch.py \
		-i helen_hdf5/prediction.hdf \
		-o /data/ \
		-p shasta_mp_helen_assembly \
		-t $(CPU)
	mv data/$(ID)/shasta_mp_helen_assembly.fa data/$(ID)/$(ID).assembly.fa
