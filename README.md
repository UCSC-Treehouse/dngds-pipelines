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
Clone this repo, download chromosome 11 sample and run all pipelines
```bash
git clone https://github.com/ucsc-upd/pipelines.git
cd pipelines
make
```
NOTE: The data directory can be a symbolic link (i.e. to a scratch location)

Several vcf files will be generated as well as various intermediate results.

## Additional Samples
To process additional samples place their fastq in data/<id>/<id>.fq.gz and call make for any specific target. For example:
```
make data/<id>/<id>.sniffles.vcf
```
