# RNA-seq Analysis Pipeline
## Datasets
- Dataset PRJNA1014743 (Jurkat T-ALL cells ± homoharringtonine, HHT)
- Blood. 2024 /  PMID: 38968151
- Homoharringtonine inhibits the NOTCH/MYC pathway and exhibits antitumor effects in T-cell acute lymphoblastic leukemia

## Folder layout
~~~
# bash
project_PRJNA1014743/
├─ raw/              # FASTQs
├─ ref/              # reference (Salmon index, GTF/FA)
├─ qc/               # FastQC & MultiQC
├─ quant/            # Salmon outputs per sample
├─ meta/             # metadata.csv, tx2gene.csv
├─ r/                # R scripts
└─ results/          # DE tables, plots, GSEA
~~~
## 0) Get the SRR runs & metadata (SRA → FASTQ)
<https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1014743>

On the PRJNA1014743 page, click “Send to” → “Run Selector” → “Run Selector” → "Download Metadata/Accession List". 
Save as **metadata.csv**
Edit a minimal metadata.csv with columns:

~~~
! [metadata.csv] ()
# CSV
sample_id,SRR,condition,replicate
CTRL_1,SRR26030909,control,1
CTRL_2,SRR26030910,control,2
HHT10_1,SRR26030907,HHT10,1
HHT10_2,SRR26030908,HHT10,2
HHT15_1,SRR26030905,HHT15,1
HHT15_2,SRR26030906,HHT15,2
~~~
## 1) Download FASTQs
**intsall fasterq-dump** 

*fasterq-dump is a tool from SRA-tools (the NCBI Sequence Read Archive toolkit)*

*fasterq-dump does not support --gzip* 

~~~
# bash
# Create a dedicated environment
conda create -n sra -c bioconda -c conda-forge sra-tools pigz
conda activate sra

# Verify installation
which fasterq-dump
fasterq-dump --version  #fasterq-dump : 3.2.1  
~~~
**Downloads the sequencing run from NCBI SRA** 
~~~
# bash
pwd # project_PRJNA1014743

cut -d, -f2 meta/metadata.csv | tail -n +2 | tr -d '\r' | while read -r SRR; do
  echo "Processing $SRR ..."
  fasterq-dump "$SRR" -e 8 --split-files -t tmp -O raw/
  pigz -p 8 raw/"${SRR}"_*.fastq   # or: gzip raw/"${SRR}"_*.fastq
done
~~~

This will produce files like:

- raw/SRR26030905/SRR26030905.sra
- raw/SRR26030906/SRR26030906.sra
- raw/SRR26030907/SRR26030907.sra
- raw/SRR26030908/SRR26030908.sra
- raw/SRR26030909/SRR26030909.sra
- raw/SRR26030910/SRR26030910.sra


~~~
find raw -type f -name "*.sra" -print0 | xargs -0 -P 6 -I{} bash -lc '
  set -e
  f="{}"; base=$(basename "$f" .sra)
  echo "[RUN] $base -> fastq_out/"

  /research/groups/yanggrp/home/glin/miniconda3/envs/sra/bin/fasterq-dump \
    --split-3 --skip-technical --threads 8 -O fastq_out "$f"

  # compress outputs (works for single-end or paired)
  shopt -s nullglob
  if command -v pigz >/dev/null; then
    pigz -p 8 fastq_out/${base}*.fastq
  else
    gzip -f fastq_out/${base}*.fastq
  fi
'
~~~
This will produce files like:

- fastq_out/SRR26030905_1.fastq
- fastq_out/SRR26030905_2.fastq
- fastq_out/SRR26030906_1.fastq
- fastq_out/SRR26030906_2.fastq
- fastq_out/SRR26030907_1.fastq
- fastq_out/SRR26030907_2.fastq
- fastq_out/SRR26030908_1.fastq
- fastq_out/SRR26030908_2.fastq
- fastq_out/SRR26030909_1.fastq
- fastq_out/SRR26030909_2.fastq
- fastq_out/SRR26030910_1.fastq
- fastq_out/SRR26030910_2.fastq

## 2) Download FASTQsDifferential Expression Design
**Since you have 3 groups**: control (PBS); HHT10; HHT15
Your DESeq2 design should be:
~~~
# r 
dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~ condition)
~~~

















