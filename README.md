## 3) Quantification 
### 3.1) Build a decoy-aware index
~~~
# bash
# ref/
mkdir -p ref

# Download (use HTTPS instead of FTP to avoid firewall issues)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

gunzip *.gz
~~~

This will produce files like:
- gencode.v44.annotation.gtf  
- gencode.v44.transcripts.fa  
- GRCh38.primary_assembly.genome.fa

~~~
# bash
# Create decoys list (chromosome headers from the genome FASTA)
grep "^>" GRCh38.primary_assembly.genome.fa | cut -d " " -f1 | sed 's/>//g' > decoys.txt
~~~
This will produce files like:
- decoys.txt 
 
~~~
# bash
# Make gentrome (transcripts + genome)
cat gencode.v44.transcripts.fa GRCh38.primary_assembly.genome.fa > gencode.v44.gentrome.fa
~~~
This will produce files like:
- gencode.v44.gentrome.fa 

~~~
conda activate sra
conda install -c bioconda salmon=1.10.3
salmon --version  #version : 1.10.3
~~~

~~~
# Build Salmon index (decoy-aware)
salmon index \
  -t gencode.v44.gentrome.fa \
  -d decoys.txt \
  -i salmon_gencode_v44_decoy \
  --gencode \
  -p 8
~~~

This will produce files like:
- Building perfect hash
- Index built successfully

### 3.2) Build a Transcript-only index and Quantify
#### 3.2.1) Build a Transcript-only index (fastest & smallest)
~~~
# In project_PRJNA1014743/ref
# Build a small, fast index (no decoys)
salmon index \
  -t gencode.v44.transcripts.fa \
  -i salmon_gencode_v44_txonly_idx \
  -p 16    # increase if your node can handle it
~~~
- Pros: tiny, fast, low RAM.
- Cons: slightly less protection against spurious mappings to unannotated genomic sequence (usually fine for bulk RNA-seq).

#### 3.2.2) Quantify
~~~
# Paired-end example
salmon quant \
  -i ref/salmon_gencode_v44_txonly_idx \
  -l A \
  -1 raw_fastq/SRR26030905_1.fastq \
  -2 raw_fastq/SRR26030905_2.fastq \
  -p 16 --validateMappings \
  -o quant/SRR26030905
~~~

### 3.3) Build a STAR index and Quantify
#### 3.3.1) Build a STAR index
First, prepare the data:
- Genome FASTA (e.g., GRCh38.primary_assembly.genome.fa)
- GENCODE annotation GTF (e.g., gencode.v44.annotation.gtf)

star_genome.sh
~~~
#!/bin/bash
#BSUB -J star_genome
#BSUB -q standard
#BSUB -n 32
#BSUB -W 24:00
#BSUB -o star_genome.%J.out
#BSUB -e star_genome.%J.err
#BSUB -R "span[hosts=1]"
# ~32 * 5.6GB ≈ 180GB total
#BSUB -R "rusage[mem=5600]"
# Hard cap in MB (>= requested total)
#BSUB -M 180000

set -euo pipefail

# initialize conda in batch shells
eval "$(conda shell.bash hook)"
conda activate sra

cd ***Sep/project_PRJNA1014743/ref
rm -rf STAR_index_gencodev44 _STARtmp
mkdir -p STAR_index_gencodev44

STAR \
  --runThreadN 32 \
  --runMode genomeGenerate \
  --genomeDir STAR_index_gencodev44 \
  --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile gencode.v44.annotation.gtf \
  --sjdbOverhang 149 \
  --limitGenomeGenerateRAM 170000000000
~~~

star_genome_build.lsf
~~~
# bash
bsub < star_genome.sh
# check the running process
bjobs
~~~

monitor 
~~~
bjobs -w -u $USER            # see when it’s RUN and on which node
bpeek -f <JOBID>             # live stdout once RUN (replace with actual ID)
tail -f star_genome.<JOBID>.out
tail -f star_genome.<JOBID>.err
tail -f STAR_index_gencodev44/Log.out
~~~
