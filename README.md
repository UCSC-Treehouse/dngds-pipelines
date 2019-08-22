# UPD Pipelines

Tooling to run the primary, secondary and tertiary pipelines for the UCSC Undiagnosed Pediatric Disease Center

## Install dependencies
To run this program, please install `Make`:
```bash
sudo apt-get update
sudo apt-get install make
```

and `Docker`:
```bash
sudo apt install docker.io
sudo systemctl start docker
sudo systemctl enable docker
```

To run `Docker` without root privilege:
```bash
sudo groupadd docker
sudo gpasswd -a $USER docker
# now log out and log in once or run:
newgrp docker
```

TIPS:
To clean all docker images to re-download  them all:
```bash
docker rmi $(docker images -a -q)
```

## Quick Start

Clone this repo, create a directory for data and run a test sample:

```bash
git clone https://github.com/ucsc-upd/pipelines.git
mkdir data
cd data
make -f ../Makefile primary secondary
```

This will generate a polished genome assembly of the input inside `GM24385.chr20_shasta_mp_helen_assembly.fa`


## Quality assessment
```
NOTE:  Shasta is non-deterministic because of dynamic load balancing.   
Because of that, repeating the same run twice under identical conditions 
does not necessarily give bit-identical results.
```

If you ran this pipeline for `GM24385.chr20` then you can assess the quality of the polished and unpolished  assembly using `pomoxis`.

Install Pomoxis:
```bash
sudo apt-get install virtualenv gcc  g++ zlib1g-dev libncurses5-dev \
python3-all-dev libhdf5-dev libatlas-base-dev libopenblas-base \
libopenblas-dev libbz2-dev liblzma-dev libffi-dev make \
python-virtualenv cmake wget bzip2

git clone --recursive https://github.com/nanoporetech/pomoxis
cd pomoxis
make install
. ./venv/bin/activate
```

Download truth assembly:
```bash
wget https://storage.googleapis.com/kishwar-helen/truth_assemblies/HG002/HG002_GRCh38_h1_chr20.fa
wget https://storage.googleapis.com/kishwar-helen/truth_assemblies/HG002/HG002_h1_chr20.bed
```

Assess the Shasta assembly:
```bash
mkdir shasta_assessment

assess_assembly \
-i shasta.fasta \
-r HG002_GRCh38_h1_chr20.fa \
-p shasta_assessment/hg002_chr20_shasta \
-b HG002_h1_chr20.bed \
-t <number_of_threads> \
-T
```

Expected output:
```bash
#  Percentage Errors
  name     mean     q10      q50      q90
 err_ont  1.181%   0.798%   1.092%   1.640%
 err_bal  1.083%   0.717%   1.008%   1.474%
    iden  0.062%   0.021%   0.056%   0.119%
     del  0.914%   0.614%   0.852%   1.212%
     ins  0.107%   0.031%   0.069%   0.188%

#  Q Scores
  name     mean      q10      q50      q90
 err_ont  19.28    20.98    19.62    17.85
 err_bal  19.65    21.44    19.97    18.31
    iden  32.04    36.74    32.49    29.26
     del  20.39    22.12    20.70    19.16
     ins  29.69    35.09    31.61    27.25
```

Then assess the polished assembly:
```bash
mkdir shasta_helen_assessment

assess_assembly \
-i GM24385.chr20_shasta_mp_helen_assembly.fa \
-r HG002_GRCh38_h1_chr20.fa \
-p shasta_helen_assessment/hg002_chr20_helen \
-b HG002_h1_chr20.bed \
-t <number_of_threads> \
-T
```

Expected output:
```bash
#  Percentage Errors
  name     mean     q10      q50      q90
 err_ont  0.526%   0.208%   0.358%   0.738%
 err_bal  0.482%   0.185%   0.333%   0.641%
    iden  0.093%   0.001%   0.033%   0.084%
     del  0.264%   0.115%   0.205%   0.384%
     ins  0.125%   0.012%   0.078%   0.210%

#  Q Scores
  name     mean      q10      q50      q90
 err_ont  22.79    26.82    24.46    21.32
 err_bal  23.17    27.33    24.78    21.93
    iden  30.32    50.00    34.80    30.75
     del  25.78    29.39    26.87    24.15
     ins  29.02    39.34    31.07    26.78
```

The assessment confirms error rate reduction of Shasta assembly from `1.181%` to `0.526%`.
