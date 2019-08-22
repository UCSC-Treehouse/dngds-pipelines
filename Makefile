ID ?= GM24385.chr20
CPU ?= 64
# reserve some CPUs while running call_consensus
HELEN_CALL_CONSENSUS_CPU ?= 56
BASE ?= $(dir $(realpath $(firstword $(MAKEFILE_LIST))))

primary:
    # this downloads the sample and converts a FASTQ file to FASTA file.
    # I think there should be a check to see if the file is FASTA or FASTQ then convert if needed.
	echo "Downloading GM24385 test sample..."
	wget -N https://lc2019.s3-us-west-2.amazonaws.com/sample_data/GM24385/GM24385.chr20.fq
	md5sum -c $(BASE)GM24385.chr20.fq.md5
	sed -n '1~4s/^@/>/p;2~4p' GM24385.chr20.fq > GM24385.chr20.fasta

secondary: shasta minimap2 samtools marginpolish helen
	
tertiary: freebayes

shasta:
    # generate Shasta assembly
	rm -rf ShastaRun
	docker run -it --rm --user=`id -u`:`id -g` --cpus="$(CPU)" -v `pwd`:/data \
		tpesout/shasta@sha256:048f180184cfce647a491f26822f633be5de4d033f894ce7bc01e8225e846236 \
		--input $(ID).fasta
	mv ShastaRun/Assembly.fasta shasta.fasta

minimap2:
    # Map the reads back to the Shasta assembly
	docker run -it --rm --user=`id -u`:`id -g` --cpus="$(CPU)" -v `pwd`:/data \
		tpesout/minimap2@sha256:5df3218ae2afebfc06189daf7433f1ade15d7cf77d23e7351f210a098eb57858 \
		-ax map-ont -t $(CPU) shasta.fasta $(ID).fq

samtools:
    # convert the mapping to to BAM
	docker run -it --rm --user=`id -u`:`id -g` --cpus="$(CPU)" -v `pwd`:/data \
		tpesout/samtools_sort:latest \
		/data/minimap2.sam -@ $(CPU)
	docker run -it --rm --user=`id -u`:`id -g` --cpus="$(CPU)" -v `pwd`:/data \
		tpesout/samtools_view@sha256:11faa9b074b3ec96f50f62133bd19f819bd5bf6ad879d913ac45955f95dd91fb \
		-hb -F 0x104 /data/samtools_sort.bam
	docker run -it --rm --user=`id -u`:`id -g` --cpus="$(CPU)" -v `pwd`:/data \
		quay.io/ucsc_cgl/samtools:1.8--cba1ddbca3e1ab94813b58e40e56ab87a59f1997 \
		index -@ $(CPU) /data/samtools_sort.bam

marginpolish:
    # generate marginPolish outputs
	rm -rf marginPolish/
	mkdir -p marginPolish
	docker run -it --rm --user=`id -u`:`id -g` --cpus="$(CPU)" -v `pwd`:/data \
		tpesout/margin_polish@sha256:de10c726bcc6af2f58cbb35af32ed0f0d95a3dc5f64f66dcc4eecbeb36f98b65 \
		/data/samtools_sort.bam /data/shasta.fasta \
		/opt/MarginPolish/params/allParams.np.human.guppy-ff-235.json -t $(CPU) -o /data/marginPolish/output -f
	mv marginPolish/output.fa marginPolish.fa

helen:
    # download the HELEN model and run call_consensus
	wget -N https://storage.googleapis.com/kishwar-helen/helen_trained_models/v0.0.1/r941_flip235_v001.pkl
	# delete previous run data (if there is any)
	rm -rf $(ID)_helen_hdf5/
	docker run -it --rm --user=`id -u`:`id -g` --cpus="$(CPU)" -v `pwd`:/data \
		kishwars/helen:0.0.1.cpu \
		call_consensus.py \
		-i /data/marginPolish/ \
		-m r941_flip235_v001.pkl \
		-o $(ID)_helen_hdf5/ \
		-p $(ID)_prediction \
		-w 0 \
		-t $(HELEN_CALL_CONSENSUS_CPU)
	docker run -it --rm --user=`id -u`:`id -g` --cpus="$(CPU)" -v `pwd`:/data \
		kishwars/helen:0.0.1.cpu \
		stitch.py \
		-i $(ID)_helen_hdf5/$(ID)_prediction.hdf \
		-o /data/ \
		-p $(ID)_shasta_mp_helen_assembly \
		-t $(CPU)

freebayes:
	echo "Downloading references..."
	wget -N http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/ENSEMBL/homo_sapiens/ENSEMBL.homo_sapiens.release-75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
	docker run -it --rm --user=`id -u`:`id -g` --cpus="$(CPU)" -v `pwd`:/data \
    quay.io/ucsc_cgl/freebayes \
			--fasta-reference /data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
			--vcf /data/freebayes.vcf \
			/data/samtools_sort.bam
