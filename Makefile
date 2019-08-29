# Sample ID - will be used to prefix all files generated
ID ?= test

# Max number of CPUs/Cores to use
CPU ?= 64

# reserve some CPUs while running call_consensus
HELEN_CALL_CONSENSUS_CPU ?= 56

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
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
		faidx /references/hg38.fa

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
# Map sample to hg38 and call variants
#

data/$(ID)/$(ID).minimap2_hg38.sam: data/references/hg38.fa data/$(ID)/$(ID).fq data/$(ID)/$(ID).fq.md5
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

data/$(ID)/$(ID).hg38.vcf: data/$(ID)/$(ID).minimap2_hg38_sorted.bam data/references/hg38.fa
	echo "Calling variants against hg38..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		--user=`id -u`:`id -g` \
		quay.io/ucsc_cgl/freebayes@sha256:b467edda4f92f22f0dc21e54e69e18bfd6dcc5cbe3292e108429e2d86034e6e5 \
		--fasta-reference /references/hg38.fa \
		--vcf /data/$(ID).hg38.vcf \
		/data/$(ID).minimap2_hg38_sorted.bam 

# http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html
data/$(ID)/$(ID).hg38_lite.vcf: data/$(ID)/$(ID).minimap2_hg38_sorted.bam data/references/hg38.fa
	echo "Calling variants on a small region against hg38..."
	docker run -it --rm --cpus="$(CPU)" -v `realpath data/$(ID)`:/data \
		-v `realpath data/references`:/references \
		--user=`id -u`:`id -g` \
		quay.io/ucsc_cgl/freebayes@sha256:b467edda4f92f22f0dc21e54e69e18bfd6dcc5cbe3292e108429e2d86034e6e5 \
		--region chr20:1000000-1010000 \
		--fasta-reference /references/hg38.fa \
		--vcf /data/$(ID).hg38_lite.vcf \
		/data/$(ID).minimap2_hg38_sorted.bam 

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
