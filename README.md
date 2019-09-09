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
Clone this repo, create samples and references directories, download references, a chromosome 11 sample and run the sniffles variant caller:
```bash
git clone https://github.com/ucsc-upd/pipelines.git
cd pipelines
mkdir -p samples references
make samples/na12878-chr11/na12878-chr11.sniffles.vcf
```
NOTE: The samples and references directories can be a symbolic links (i.e. to a scratch location or into a shared file system)

This will take approximately 72 minutes using 32 cores and genereate the following output in samples/na12878-chr11:
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
To process additional samples place their fastq in samples/<id>/<id>.fq.gz and call make for any specific target. For example:
```
make samples/<id>/<id>.sniffles.vcf
```

## Other Targets
```
make samples/na12878-chr11/na12878-chr11.svim.vcf
make samples/na12878-chr11/na12878-chr11.clairvoyant.vcf
```
