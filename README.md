# RatGTEx Pipeline

This is the code used to process data for the [RatGTEx Portal](https://ratgtex.org). It is built on [Snakemake](https://snakemake.github.io/), a Python-based framework for reproducible data analysis. Since things might not work the first time with a new dataset, it makes it easy to process data iteratively. You can run part of the pipeline on a subset of the data, and then run it on the full dataset without regenerating existing files.

The main steps of the pipeline are:
0. Download RNA-Seq FASTQ files from the [SRA](https://www.ncbi.nlm.nih.gov/sra) and genotypes from RatGTEx, and phase the genotypes with [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html).
1. Align RNA-Seq reads using [STAR](https://github.com/alexdobin/STAR).
2. Checking for and fixing sample mixups.
3. Quantify RNA phenotypes using [Pantry](https://github.com/PejLab/Pantry).
4. Map cis-eQTLs and trans-eQTLs using [tensorQTL](https://github.com/broadinstitute/tensorqtl) in various modes. Some of these are built into the Pantry Pheast module, while others are run in the main RatGTEx snakemake pipeline.
5. Compute TWAS weights for all RNA phenotypes using the Pantry Pheast module. These are used for the Rat TWAS Hub.
6. Calculate cis-eQTL effect size (allelic fold change) using [aFC.py](https://github.com/secastel/aFC).

## Setup

### Conda environment (for launching Snakemake only)

Install a conda-like package manager (I recommend [miniforge](https://github.com/conda-forge/miniforge)) and add the bioconda channel.

Create the ratgtex environment, which only needs to run Snakemake itself:

```shell
conda env create -n ratgtex --file environment.yml
conda activate ratgtex
```

All of the actual bioinformatics tools run inside containers (below), not from this environment.

### Containers

The pipeline runs its tools inside [Apptainer](https://apptainer.org/) (formerly Singularity) images rather than a conda environment, which avoids conda pulling binaries that don't run on an older cluster's glibc. The container definitions are in `containers/`:

- `sratools.def` — the SRA Toolkit (`prefetch`, `fasterq-dump`) for downloading FASTQs.
- `bioinfo.def` — general tools for alignment, QC, genotype processing, and phasing: STAR, samtools, bcftools, htslib (`tabix`/`bgzip`), plink2, plink, gatk4, Beagle, and Python with pandas/numpy/scipy/pysam/pyyaml/fastparquet/bx-python.
- `tensorqtl.def` — a GPU (CUDA/PyTorch) image for the tensorQTL rules. Only needed when the QTL steps are re-enabled; build it on a machine matching the target CUDA version.

Build the images (this is also done in `scripts/setup/setup_v4.sh`):

```shell
mkdir -p images
apptainer build --ignore-subuid --ignore-fakeroot-command images/sratools.sif containers/sratools.def
apptainer build --ignore-subuid --ignore-fakeroot-command images/bioinfo.sif containers/bioinfo.def
# apptainer build --ignore-subuid --ignore-fakeroot-command images/tensorqtl.sif containers/tensorqtl.def
```

Each rule has a `container:` directive pointing at one of `images/*.sif`, so run Snakemake with `--use-singularity` (as in `run_pipeline.sh`). Do not pass `--use-conda` — the two aren't combined here. Bind any paths the jobs need to read/write with `--singularity-args "-B /your/path"`.

### Other software

To get `aFC.py`:

```
cd tools
git clone git@github.com:secastel/aFC.git
```

### Snakemake profile

When you run snakemake, you specify a profile that determines how steps get run. Here is the config file I use on a computing cluster with slurm scheduling:

`~/.config/snakemake/slurm/config.yaml`:

```yaml
executor: slurm
default-resources:
  runtime: "4h"
  mem_mb: 8000
  slurm_partition: "our-cpu-partition"
  slurm_account: "our-cpu-account"
set-resources:
  tensorqtl_cis:
    slurm_partition: "our-gpu-partition"
    slurm_account: "our-gpu-account"
    slurm_extra: "'--gres=gpu:1'"
  tensorqtl_cis_independent:
    slurm_partition: "our-gpu-partition"
    slurm_account: "our-gpu-account"
    slurm_extra: "'--gres=gpu:1'"
  tensorqtl_nominal:
    slurm_partition: "our-gpu-partition"
    slurm_account: "our-gpu-account"
    slurm_extra: "'--gres=gpu:1'"
  tensorqtl_trans:
    slurm_partition: "our-gpu-partition"
    slurm_account: "our-gpu-account"
    slurm_extra: "'--gres=gpu:1'"
resources:
  sra_downloads: 3
use-singularity: true
latency-wait: 60
```

This uses snakemake v8 or higher and the `snakemake-executor-plugin-slurm` plugin. The `tensorqtl` steps should be run on GPU for reasonable runtime, so they are specified by name here to override the default resources. Adjust as needed for your cluster. Additional resources are specified within some of the snakemake rules, which are passed to slurm when those jobs are submitted. Alternatively, you can run snakemake on an interactive node.

The `resources: sra_downloads: 3` line caps concurrent SRA downloads at 3. The `sra_fastq_paired` and `sra_fastq_single` rules each consume `sra_downloads=1`, and this sets the total pool. Change the number to allow more or fewer simultaneous downloads (or override at run time with `--resources sra_downloads=N`). If you run Snakemake without a profile, pass it on the command line, e.g. `snakemake --resources sra_downloads=3 ...`.

#### `config.yaml`

This configuration file contains parameters about the datasets and reference data and is used by snakemake. Parameters for each tissue are grouped under tissue names, and the `run` list controls which tissues are processed, e.g.:

```yaml
version: "v4"
ref_genome: "ref/GCF_036323735.1_GRCr8_genomic.chr.fa"
ref_anno: "ref/GCF_036323735.1_GRCr8_genomic.chr.gtf"

run:
  - IL
  - LHb

tissues:
  # IL, LHb, NAcc1, OFC, and PL1 datsets
  IL:
    read_length: 100
    fastq_path: "fastq/IL_LHb_NAcc_OFC_PL"
    paired_end: false
    geno_dataset: "ratgtex_v4_round11_2"
  LHb:
    read_length: 100
    fastq_path: "fastq/IL_LHb_NAcc_OFC_PL"
    paired_end: false
    geno_dataset: "ratgtex_v4_round11_2"
...
```

Some of these are inherent to the data, while others, e.g. `fastq_path`, may need to be edited to point to the correct location.

### Dataset-specific input files

#### FASTQ files (downloaded from the SRA)

The RNA-seq FASTQ files are downloaded from the SRA by the pipeline. RatGTEx
sample IDs are mapped to SRA accessions in the sample table at
`https://ratgtex.org/data/{version}/ref/RatGTEx_samples.{version}.tsv`, whose
columns are `rat_id`, `tissue`, `GEO_accession`, `BioSample_accession`, and
`SRA_accession`.

Download that table and run `scripts/setup/setup_sra.py` to generate the
per-tissue input files:

```shell
wget https://ratgtex.org/data/v4/ref/RatGTEx_samples.v4.tsv
python3 scripts/setup/setup_sra.py \
    --samples RatGTEx_samples.v4.tsv \
    --version v4 \
    --tissue Liver
```

This resolves each `SRA_accession` (usually an SRX experiment) to its run
accession(s) and library layout via the public ENA portal API, then writes
`rat_ids.txt`, `fastq_map.txt`, and a `sra_runs.tsv` record for the tissue. When
Snakemake runs, the `sra_fastq_paired` / `sra_fastq_single` rules download each
run with `prefetch` + `fasterq-dump` into `fastq/{tissue}/`. FASTQ files can be
single-read or paired-end, and there can be more than one run per sample; all
will be aligned into one BAM file.

The Liver dataset is paired-end for every sample. Other datasets may contain
single-end (or a mix of single- and paired-end) samples; `setup_sra.py` records
the correct layout per run, and the pipeline handles both.

#### `{version}/{tissue}/fastq_map.txt`

A tab-delimited file with no header containing the paths to each FASTQ file and the rat IDs they correspond to. Or, for paired-end reads, each row contains the first FASTQ path, second FASTQ path, and rat ID per file pair. It is generated by `setup_sra.py`, but the format is documented here for reference.
- If multiple files map to the same ID, i.e. the ID appears in multiple rows, reads from those files will be aligned into one BAM file.
- Paths are relative to the `fastq_path` parameter in `config.yaml` (e.g. `fastq/Liver`), which is where the SRA download rules write the FASTQs.
- Any listed files whose rat IDs are not in `{version}/{tissue}/rat_ids.txt` will be ignored.
- If a tissue includes a mix of single-end and paired-end samples, omit `paired_end` from the tissue config so the pipeline infers pairing per row.

#### `{version}/{tissue}/rat_ids.txt`

A file listing the rat IDs for the dataset, one per line. This list determines which samples are included in the processing. It is generated by `setup_sra.py` (only rats with at least one resolvable SRA run are included).

#### `geno/{dataset}.vcf.gz`

A VCF file containing the unphased genotypes for one or more tissues. If multiple tissues came from the same project and have overlapping sets of individuals, they use the same VCF file. The pipeline downloads the all-rat genotype VCF from the URL given in `config.yaml` under `geno_url` (the `download_genotypes` rule), then the `process_genotypes` rule ensures REF alleles match the reference genome and keeps only biallelic SNPs — the same normalization documented in `scripts/setup/genotypes_{version}.sh`. Specify the dataset name in the pipeline config file as described below.

#### `geno/{dataset}.phased.vcf.gz`

Phased genotypes produced by the `phase_genotypes` rule using Beagle
(reference-free / population phasing, since there is no rat haplotype reference
panel). An optional genetic map per dataset can be supplied in `config.yaml`
under `geno_map`; without one Beagle assumes a uniform 1 cM/Mb rate. These
phased genotypes are built by default (they are listed in the `all` rule) for
use in downstream analyses.

## Running

Edit `config.yaml` in this directory so that the tissue(s) you want to process are present in the `tissues` section with correct parameters, and are listed in the `run` section. Unlike the Snakemake config file, which specifies how jobs are run, this one contains parameters for the tissues/datasets such as read length and directory where FASTQ files can be found.

### QC

#### Pre-run checks

Before running Snakemake, run `python3 scripts/setup/init_check.py {version} {tissue}`, which checks the input data and config for issues.

#### Sample mixup checks

The way to do sample mixup testing is to generate the mixup checking outputs using Snakemake, which will generate the BAM files as dependencies if needed. Examine the outputs to identify samples that need to be relabeled (e.g. if two labels get swapped) or removed.

- To relabel a sample, edit the rat ID in the last column of `fastq_map.txt` for all of its FASTQ files so that its BAM file gets labeled correctly. You'll then need to regenerate the BAM file since it will now use the correct VCF individual as input to STAR.
- To remove a sample, remove its ID from `rat_ids.txt` and delete its BAM and any other generated files.

Before removing samples, run the second stage of sample mixup checking, which tests the RNA-seq samples that still don't have matches against 6000+ rat genotypes to see if a match can be found. To do this, list the mismatched samples in `{version}/{tissue}/qc/samples_without_matches.txt`, along with an OK sample as a positive control (if that sample is included in the all-rat VCF). Then generate `{version}/{tissue}/qc/all_rats_summary.tsv` and use any additional matches found. This will probably require adding the new matching genotypes to the VCF file (see `scripts/setup/genotypes_{version}.sh`).

### Continue

After these corrections, continue with the pipeline.

You may want to run a subset of the heavy raw data processing steps first, then move on once those are done. E.g. add the first 10 BAM files to the first rule (called 'all') in `Snakefile` and generate them:

`snakemake --profile slurm -j10`

Use the `-n` dry run tag to make sure things seem to be set up correctly before running.

### Merging same-tissue datasets

Sets of samples from the same tissue, collected by different studies, can be run individually and then run as a merged dataset. These merges (e.g. NAcc1 + NAcc2 + NAcc3 = NAcc) are specified in the config. Keep the merged tissue names commented out in the `run` list at first to run the individual datasets, at least before running Pantry phenotyping. Run `scripts/setup/setup_merged_tissues.sh` to set up the merged directories. Then, uncomment the merged tissue names in the `run` list and run them.

## Help

There are likely a number of issues remaining with this pipeline, so email me or file an issue on this GitHub repository if anything isn't working, and I'll be happy to fix it.
