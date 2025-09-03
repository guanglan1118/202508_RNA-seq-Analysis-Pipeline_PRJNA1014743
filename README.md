# RNA-seq Analysis Pipeline
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
On the PRJNA1014743 page, click “SRA Experiments (6)” → “Run Selector” → Download the table (CSV). Save as meta/sra_run_table.csv.
Edit a minimal metadata.csv with columns:

*CSV*
~~~
sample_id,condition,replicate,SRR
CTRL_1,control,1,SRRXXXXXXXX
CTRL_2,control,2,SRRXXXXXXXX
CTRL_3,control,3,SRRXXXXXXXX
HHT_1,HHT,1,SRRXXXXXXXX
HHT_2,HHT,2,SRRXXXXXXXX
HHT_3,HHT,3,SRRXXXXXXXX
~~~
