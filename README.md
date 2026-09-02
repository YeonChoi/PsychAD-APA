# PsychAD-APA

Alternative polyadenylation (APA) analysis of single-nucleus RNA-seq from human
prefrontal cortex, as used in the PsychAD study.

The pipeline calls and quantifies polyadenylation sites (PAS) per cell type and
donor with [SCAPTURE](https://github.com/YangLab/SCAPTURE), then tests them for
differential expression (dreamlet) and differential usage (crumblr).

 **Scope of this repository.** This is a cleaned version of the analysis.
 The exact code and parameters used to produce the published figures
 are archived separately at [Zenodo DOI, TBD].

---

## Quick start

Install SCAPTURE, subread, samtools separately — see [Prerequisites](#prerequisites).

```bash
# 1. Python/R environment
conda env create -f conda_base_APA.no-builds.yaml
conda activate base_APA

# 2. Edit config.yaml: gene_h5ad, barcode_table, scapture_env
#    (scapture_env = name of the conda env you installed SCAPTURE into)

# 3. Quantify PAS
cd 2_run_SCAPTURE_quant && snakemake -j --use-conda
```

---

## Repository layout

Not everything in this repository is ready-to-run, because everyone has different data structures / questions.
The main pipeline is in `2_run_SCAPTURE_quant`, which generates a single-cell PAS `.h5ad` file and pseudobulk `.csv.gz` files.

🟢 use as-is &nbsp;&nbsp; 🟡 example - read and adapt to your own data

| Component | Description | Ready-to-use |
|---|---|:---:|
| `config.yaml` | All input/output paths | 🟡 |
| `conda_base_APA.no-builds.yaml` | Python/R environment for Snakemake and the notebooks | 🟢 |
| `helpers/` | Shared Python utilities imported by the notebooks and Snakefiles | 🟢 |
| `ref/` | PAS annotation called on PsychAD (BED12 + per-PAS metadata) | 🟢 |
| `1_barcode_tables/` | **Stage 1** - build the cell-barcode / BAM-path table from a gene-level h5ad | 🟡 |
| `2_run_SCAPTURE_quant/` | **Stage 2 (main)** - SCAPTURE PASquant per Channel (Snakemake), then merge to a PAS h5ad and pseudobulk (notebooks) | 🟢 |
| `3_run_dreamlet/` | **Stage 3** - differential expression with dreamlet (R Markdown) | 🟡 |
| `4_run_crumblr/` | **Stage 4** - differential PAS usage on within-gene CLR ratios with crumblr (R Markdown) | 🟡 |
| `5_calc_dWUI/` | **Stage 5** - combine DE and DU results into a per-gene weighted usage index shift | 🟡 |
| `99_run_SCAPTURE_annotation/` | How `ref/` was produced: PAS calling, DB/DeepPASS filtering, cross-gene merging, renaming | 🟡 |

---


## Prerequisites

### Must be on `$PATH`

Install SCAPTURE in a conda environment, following the
[SCAPTURE instructions](https://github.com/YangLab/SCAPTURE). Then set
`scapture_env` in `config.yaml` to that environment's name (default: `SCAPTURE_env`).

| tool | version used here |
|---|---|
| scapture | v1.1 (optionally + `modifications_on_scapture.diff`, see below) |
| bedtools | 2.26.0 (pinned by SCAPTURE) |
| featureCounts (subread) | 1.6.3 |
| samtools | 1.21 |
| umi_tools | 1.1.0 |

featureCounts is not bundled by SCAPTURE's environment; add it to `$PATH`
yourself. Its version affects read assignment (`-M -O --largestOverlap -s 1`),
so we recommend pinning 1.6.3.

SCAPTURE is GPLv3 and is **not** redistributed here.

`modifications_on_scapture.diff` records the small changes we made to SCAPTURE for
our runs: thread counts, a looser PAS overlap cutoff, and a mapping-quality
filter in the UMI counting step. **Applying it is optional** — the pipeline runs
on stock SCAPTURE, and only the last two affect results. Apply it to a clean
SCAPTURE checkout if you want to match our runs exactly.

### Python / R

```bash
conda env create -f conda_base_APA.no-builds.yaml
conda activate base_APA
```

Everything in stages 1, 2 and 5 runs in this environment. Stages 3 and 4 are
R Markdown and additionally need `dreamlet`, `crumblr`, `zellkonverter` and
`SingleCellExperiment`.

---

## Inputs

| input | what it is |
|---|---|
| gene-level `.h5ad` | one file, `.obs` supplies cell barcodes, `Channel`, `individualID` and cell-type labels (`class`/`subclass`/`subtype`) |
| aligned BAMs | one per sequencing batch, coordinate-sorted and indexed, with `CB`/`UB` tags (STARsolo or CellRanger) |
| PAS annotation | `ref/all_passed_merged_PAS-ageXclass.bed` — 108,768 PAS called on the PsychAD cohort, shipped here |

The `.obs` index is expected to look like
`Donor123-rep2_AAACCCAAGAGGGTCT`, i.e. `<Channel>_<barcode>`.

Point `config.yaml` at your own files:

```yaml
dataset_id: # For example, 'apa_test'
gene_h5ad:  # path of gene-level single-cell h5ad file. 
pas_anno_bed: ref/all_passed_merged_PAS-ageXclass.bed
pas_info_tsv: ref/PAS_merged_evaluated_with_name_ageXclass.tsv
barcode_table: 1_barcode_tables/apa_test_bam_path_and_barcodes.parquet
tools_dir: helpers/
```

Paths are relative to `config.yaml` itself (the repository root); absolute paths
also work.

---

## Running

### 1. Barcode table

The purpose of this section is to generate a table in the format below. If you
already have one, you can skip this.

| Channel | snRNA_bam | cell_barcode |
|---|---|---|
| Donor123-rep1 | pathA/to_B.bam | AAACCCAAGAGGGTCT |
| Donor123-rep1 | pathA/to_C.bam | CCCCCTAAGAAAGTCT |
| Donor456-rep1 | pathD/to_E.bam | TTTCCCGAGGGTCTGG |

The example notebook reads `.obs` from the gene `.h5ad`, derives one BAM path per
`Batch`, and writes `barcode_table` in `.parquet` format. Other tabular formats
such as `.csv` and `.tsv` are also accepted.

### 2. PAS quantification

#### Snakemake

Make sure all file paths in `config.yaml` are correct, then:

```bash
# for single machine
cd 2_run_SCAPTURE_quant && snakemake -j --use-conda
# for LSF managed computing cluster
cd 2_run_SCAPTURE_quant && ./run_snake_mount-sinai_lsf.sh
```

One SCAPTURE `PASquant` job per `Channel`. Outputs land in
`counts/{channel}/{channel}.KeepCell.UMIs.tsv.gz`, plus a per-channel
`.KeepPAS.metadata` (the PAS annotation after SCAPTURE's overlap filtering and
renaming — identical across channels).

`run_snake_mount-sinai_lsf.sh` and `lsf.yaml` are written for Mount Sinai LSF. Modify according to your environment.

**Runtime.** For reference, one run over ~400 Channels (~2.2M nuclei) took
**5.5 h wall-clock on 100 LSF nodes**, about **2,300 CPU-hours** in total.
Per Channel `scapture_pas_quant` takes ~5.6 CPU-hours on 6 threads, but varies
widely (9–175 min wall-clock).

#### ipynb

In `2_run_SCAPTURE_quant/ipynb/`:

- `1_make_h5ad_and_pseudobulk_tables.ipynb` — per-channel counts → merged single-cell PAS `.h5ad` → pseudobulk CSVs per cell type
- `2_Sum_Channel_counts_to_merged_level.ipynb` — sum Channels up to donor level
- `QC_Check_read_counts.ipynb`, `QC_Check_gene-PAS_corr.ipynb` — QC examples

#### What you get

- `counts/{channel}/{channel}.KeepCell.UMIs.tsv.gz`
  PAS x cell UMI counts, one file per Channel
- `stats/merged_count.txt`
  reads assigned / unassigned per Channel
- `singlecell_counts/*_PAS.h5ad`
  all Channels merged into one single-cell PAS object
- `singlecell_counts/PAS_info_for_rowdata_all.csv.gz`
  PAS coordinates, gene names and PAS groups (the `.var` table)
- `psb/count_no_cutoff/{class,subclass,subtype}/PAS_read_count_{celltype}_{level}.csv.gz`
  pseudobulk PAS counts, per cell type, at Channel and donor level

### 3. Differential expression (dreamlet)

> Stages 3–5 are provided as **reference implementations**. They still carry the
> directory layout of our HPC and model formulas specific to our hypotheses, and
> are meant to be read and adapted, not executed unchanged.

See the [dreamlet](https://gabrielhoffman.github.io/dreamlet/) documentation for more details.

`3_run_dreamlet/` — `run_dreamlet_voom-{Gene,PAS}.Rmd` builds pseudobulk and runs voom; the `AD_vs_Ctrl/`, `BRAAK/`, `CERAD/` subdirectories run one contrast each.

### 4. Differential usage (crumblr)

See the [crumblr](https://gabrielhoffman.github.io/crumblr/index.html) and [dream](https://gabrielhoffman.github.io/variancePartition/articles/dream.html) documentation for more details.

`4_run_crumblr/` — `makeCrumblr_pergene_PAS_ratio.Rmd` converts PAS counts to
within-gene compositional ratios (CLR); the contrast subdirectories test them.

The examples here perform only the all-PAS CLR test; for other stratifications,
such as within the 3′ UTR or per genomic region, see the full project on Zenodo.

### 5. dWUI

`5_calc_dWUI/` — combines the DE and DU results into a per-gene weighted usage
index shift.

---

### 99. PAS annotation

The two files in `ref/` were produced by this stage. You only need it if you want
to call PAS on your own cohort instead of reusing ours.

> This stage documents **how the annotation was made**, not a turnkey pipeline.
> It needs raw BAMs, a reference genome and several public polyA databases, and
> one step depends on long-read data we cannot redistribute (see below).

| Step | What it does |
|---|---|
| `1_make_pas/1_anno/` | `scapture -m annotation` — converts a GTF and genome FASTA into the exonic / intronic / 3′-extended gene models that PAScall searches |
| `1_make_pas/2_subsample_cell_types/` | Reads `.obs` from the gene h5ad and samples N cells per age X cell-class group, writing per-BAM barcode lists |
| `1_make_pas/3_pl_ageXclass/` | Splits BAMs by group, merges per group, subsamples every group to equal depth, then runs `scapture -m PAScall` |
| `2_filter_pas/` | Intersects the called peaks with four public polyA databases, then keeps a PAS if it is supported by two databases, **or** only one database and passes DeepPASS |
| `3_merge_pas/` | Merges PAS shared between overlapping genes (recursive interval merge) and reassigns each PAS back to every gene it maps to |
| `4_rename_quantified_PAS/` | Assigns the final PAS names and the gene / PAS-group mapping, generating `PAS_merged_evaluated_with_name_ageXclass.tsv` |

Stratifying by age X cell class (`pl_ageXclass`) and subsampling each group to
equal read depth keeps PAS discovery from being dominated by the most abundant
cell types. Substitute your own grouping if that is not what you need.

**Reference data required** — none of it is shipped here:

| file | source |
|---|---|
| GRCh38 primary assembly FASTA + chrom sizes | GENCODE / Ensembl |
| Ensembl 104 GTF | Ensembl |
| `SupTab_KnownPASs_fourDBs.txt` | bundled with SCAPTURE |
| GENCODE v38 `polyA_site` | GENCODE |
| PolyA-Seq, merged hg38 | published atlas |
| PolyA_DB 3 | published atlas (hg19, lifted over to hg38) |
| PolyASite 2.0 and 3.0 | polyasite.unibas.ch |

All of them are converted to Ensembl-style chromosome names (`1`, not `chr1`)
before use.

**What you cannot reproduce exactly**

- When one PAS maps to several overlapping genes, `3_merge_pas` breaks the tie
  using PsychAD pseudobulk CPM and long-read CPM. The long-read table is not
  redistributed, so that column is unavailable;
  drop the term or supply your own long-read support to run this step.
- PsychAD BAMs are controlled-access (see **Data availability**).


---

## Data availability

Processed single-nucleus data for the PsychAD cohort are available from the
Synapse portal (https://www.synapse.org) under accession `syn74739724`. This
pipeline starts from aligned BAMs, so to run stage 2 on PsychAD you will need to
generate your own BAMs from the PsychAD raw reads in the AD Knowledge Portal
(https://adknowledgeportal.org).

## Citation

Please cite the following:
```
[TBD]
```

## License

MIT — see [LICENSE](LICENSE).
