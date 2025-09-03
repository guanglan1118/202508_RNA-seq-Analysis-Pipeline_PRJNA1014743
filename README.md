me se# RNA-seq Analysis Pipeline
## Datasets
- dataset PRJNA1014743 (Jurkat T-ALL cells ± homoharringtonine, HHT)
- Blood. 2024 /  PMID: 38968151
- Homoharringtonine inhibits the NOTCH/MYC pathway and exhibits antitumor effects in T-cell acute lymphoblastic leukemia

## Folder layout
*bash*
~~~
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
On the PRJNA1014743 page, click “Send to” → “Run Selector” → “Run Selector” → "Download Metadata/Accession List". 
Save as **metadata.csv**
Edit a minimal metadata.csv with columns:

*CSV*
~~~
sample_id,SRR,condition,replicate
CTRL_1,SRR26030909,control,1
CTRL_2,SRR26030910,control,2
HHT10_1,SRR26030907,HHT10,1
HHT10_2,SRR26030908,HHT10,2
HHT15_1,SRR26030905,HHT15,1
HHT15_2,SRR26030906,HHT15,2
~~~
## 1) Download FASTQs
*bash*
~~~
mkdir -p raw/
cut -d, -f2 metadata.csv | tail -n +2 | while read SRR; do
  fasterq-dump $SRR --split-files --gzip -O raw/
done
~~~
This will produce files like (for paired-end data):

- raw/SRR26030905_1.fastq.gz

- raw/SRR26030905_2.fastq.gz

## 2) Download FASTQsDifferential Expression Design


















