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
