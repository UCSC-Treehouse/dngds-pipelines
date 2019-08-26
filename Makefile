ID ?= test
CPU ?= 64
# reserve some CPUs while running call_consensus
HELEN_CALL_CONSENSUS_CPU ?= 56
APP_PATH ?= $(dir $(realpath $(firstword $(MAKEFILE_LIST))))

#
# Download references
#

references/hg38.fa:
	# Download references any of the pipelines require
	mkdir -p references
	wget -N -P references https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
	wget -N -P references http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
	chmod +x references/twoBitToFa
	references/twoBitToFa references/hg38.2bit references/hg38.fa

references/hg38.fai:
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		-v `realpath ../references`:/references \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
			faidx /references/hg38.fa

#
# Download a test sample
#

sample.fq:
	# Download a reference quality nanopore fastq pruned to chromosome 20
	wget -N -O sample.fq https://lc2019.s3-us-west-2.amazonaws.com/sample_data/GM24385/GM24385.chr20.fq
	sed -n '1~4s/^@/>/p;2~4p' sample.fq > sample.fa


#
# Map sample to hg38 and call variants
#

minimap2.hg38.sam: ../references/hg38.fa sample.fq
	# Map the reads to HG38
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		-v `realpath ../references`:/references \
		tpesout/minimap2@sha256:5df3218ae2afebfc06189daf7433f1ade15d7cf77d23e7351f210a098eb57858 \
			-ax map-ont -t $(CPU) /references/hg38.fa sample.fq > minimap2.hg38.sam

minimap2.hg38.bam:
	# Convert mapped sam to bam
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
			view -S -b /data/minimap2.hg38.sam -o /data/minimap2.hg38.bam

minimap2.hg38.sorted.bam:
	# Convert mapped sam to bam
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
			sort /data/minimap2.hg38.bam -o /data/minimap2.hg38.sorted.bam

minimap2.hg38.sorted.bam.bai:
	# Convert mapped sam to bam
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
			index /data/minimap2.hg38.sorted.bam

freebayes.hg38.vcf:
	# Call variants against hg38
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		--user=`id -u`:`id -g` \
		-v `realpath ../references`:/references \
    quay.io/ucsc_cgl/freebayes@sha256:b467edda4f92f22f0dc21e54e69e18bfd6dcc5cbe3292e108429e2d86034e6e5 \
			--fasta-reference /references/hg38.fa \
			--vcf /data/freebayes.hg38.vcf \
			/data/minimap2.hg38.sorted.bam 

#
# De novo assemble and polish the sample
#

shasta:
	# Generate Shasta assembly
	rm -rf ShastaRun
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		tpesout/shasta@sha256:048f180184cfce647a491f26822f633be5de4d033f894ce7bc01e8225e846236 \
		--input sample.fa
	mv ShastaRun/Assembly.fasta shasta.fa

minimap2:
	# Map the reads back to the Shasta assembly
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		tpesout/minimap2@sha256:5df3218ae2afebfc06189daf7433f1ade15d7cf77d23e7351f210a098eb57858 \
			-ax map-ont -t $(CPU) shasta.fa sample.fq > minimap2.shasta.sam

samtools:
	# Convert the mapping to BAM
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
			view -S -b /data/minimap2.shasta.sam -o /data/minimap2.shasta.bam
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
			sort /data/minimap2.shasta.bam -o /data/minimap2.shasta.sorted.bam
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		quay.io/ucsc_cgl/samtools@sha256:2abed6c570ef4614fbd43673ddcdc1bbcd7318cb067ffa3d42eb50fc6ec1b95f \
			index /data/minimap2.shasta.bam

marginpolish:
	# Generate marginPolish outputs
	rm -rf marginPolish/
	mkdir -p marginPolish
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		tpesout/margin_polish@sha256:de10c726bcc6af2f58cbb35af32ed0f0d95a3dc5f64f66dcc4eecbeb36f98b65 \
		/data/samtools_sort.bam /data/shasta.fasta \
		/opt/MarginPolish/params/allParams.np.human.guppy-ff-235.json -t $(CPU) -o /data/marginPolish/output -f
	mv marginPolish/output.fa marginPolish.fa

helen:
	# Download the HELEN model and run call_consensus
	wget -N https://storage.googleapis.com/kishwar-helen/helen_trained_models/v0.0.1/r941_flip235_v001.pkl
	delete previous run data (if there is any)
	rm -rf helen_hdf5/
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		kishwars/helen@sha256:ac3504a6450c57b138a51652ebfff39cf1f5a0435ec995fcf137c7a1978cbedd \
		call_consensus.py \
		-i /data/marginPolish/ \
		-m r941_flip235_v001.pkl \
		-o helen_hdf5/ \
		-p prediction \
		-w 0 \
		-t $(HELEN_CALL_CONSENSUS_CPU)
	docker run -it --rm --cpus="$(CPU)" -v `pwd`:/data \
		kishwars/helen@sha256:ac3504a6450c57b138a51652ebfff39cf1f5a0435ec995fcf137c7a1978cbedd \
		stitch.py \
		-i helen_hdf5/prediction.hdf \
		-o /data/ \
		-p shasta_mp_helen_assembly \
		-t $(CPU)
