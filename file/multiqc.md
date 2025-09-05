# multiqc report 
## 1) General Statistics
- Dups (%) → Percentage of duplicate reads. High duplication could mean PCR bias or low library complexity. (20~40%)

- GC (%) → GC content (percentage of G and C bases). For human samples, ~50% is normal.

- Seqs → Total number of reads per sample. Yours are ~24–25 million reads each, which is good depth

## 2) FastQC
### 2.1) Sequence Counts
- Total sequencing depth: ~12–17 million reads per sample.

- Unique reads: ~40–45% per sample.

- Duplicate reads: ~55–60% per sample, consistent with the duplication table.

### 2.2) Sequence Quality Histograms 
- All samples show mean Phred scores around 33–35 at the start, with a slight decrease toward the read ends (~31–32 at 100 bp).

- The entire curve stays well within the green zone (Q ≥ 30) across all read positions.

- This means that your sequencing reads are of excellent quality across their full length.
- 
### 2.3) Per Sequence Quality Scores
- The distribution (green line) is essentially a sharp peak at the far right side of the plot.

- No visible distribution in the yellow or red zones.

- This means that nearly all your sequencing reads are of very high quality.

- There is no subset of poor-quality reads.

- This is consistent with the previous “per-base quality” plot — excellent overall sequencing.

## 3) GC Content

- Per Sequence GC Content → Shows distribution of GC content across reads.

- Ideally, this should look like a normal (bell-shaped) distribution centered around 50%.

- Your samples cluster close to 49–50%, which looks normal 


## 4) Sequence Length

- Shows how long the reads are (e.g., 100 bp, 150 bp).

- Yours are consistent across samples, which is good 


## 5) Duplication & Overrepresented Sequences

- Duplication Levels → Low duplication is better (yours are moderate, not alarming).

- Overrepresented sequences → Less than 1% in all samples, so no major contamination 


## 6) Adapter Content

- Checks for leftover adapter sequences.

- Your samples show very low adapter contamination, which is good 


## 7) Status Checks

- Each QC category (like base quality, GC content, adapter content) is marked:

- ✅ Green = normal

- 🟧 Orange = slightly unusual

-❌ Red = very unusual

- From your report, most samples are green, meaning the sequencing quality is overall good 


## 8) ✅ Bottom Line

- Your sequencing data looks high quality:

- ~25M reads per sample

- ~50% GC content (normal for human)

- Low adapter contamination

- <1% overrepresented sequences

- Quality scores are consistently high

So, your dataset is ready for downstream analysis (like alignment, quantification, or variant calling).





