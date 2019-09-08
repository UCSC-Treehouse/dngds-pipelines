# UPD Pipelines
Tooling to run primary, secondary and tertiary pipelines for the UCSC Undiagnosed Pediatric Disease Center

## Requirements
* linux
* make
* docker
* 50G disk space for the chromosome 11 sample
* 256GB memory (?)
* 32+ cores (?)
* GPU for clairvoyant

## Quick Start
Clone this repo, create data directory and download references and a chromosome 11 sample and run the sniffles variant caller:
```bash
git clone https://github.com/ucsc-upd/pipelines.git
cd pipelines
mkdir data
make data/na12878-chr11/na12878-chr11.sniffles.vcf
```
NOTE: The data directory can be a symbolic link (i.e. to a scratch location)

This will take approximately 72 minutes using 32 cores and genereate the following output in data/na12878-chr11:
```
1.3K Sep  8 11:00 minimap2.log
4.7G Sep  8 11:10 na12878-chr11.bam
3.6G Sep  8 10:25 na12878-chr11.fq.gz
  73 Sep  8 10:25 na12878-chr11.fq.gz.md5
5.8G Dec 21  2016 na12878-chr11.original.bam
 11G Sep  8 11:00 na12878-chr11.sam
835K Sep  8 11:29 na12878-chr11.sniffles.vcf
4.8G Sep  8 11:24 na12878-chr11.sorted.bam
```

## Additional Samples
To process additional samples place their fastq in data/<id>/<id>.fq.gz and call make for any specific target. For example:
```
make data/<id>/<id>.sniffles.vcf
```

## Other Targets
```
make data/na12878-chr11/na12878-chr11.svim.vcf
make data/na12878-chr11/na12878-chr11.clairvoyant.vcf
```
