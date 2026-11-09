[MonopogenLite](https://github.com/swvanderlaan/MonopogenLite)<img align="right" height="200" src=images/MonopogenLite_logo_print.png>
============
[![Languages](https://skillicons.dev/icons?i=bash,r,py)](https://skillicons.dev) 

**MonopogenLite** _Germline SNV calling and phasing from single-cell sequencing data (for macOS Sequoia and Linux Rocky8)._

This is a fork of the original [`Monopogen`](https://github.com/KChen-lab/Monopogen) which works with python (3.9+), and `samtools` and `bcftools` (v1.21), `vcftools` (v0.1.16), and `tabix` (`htslib`, v1.21), as well as in the context of `Rocky8` (Linux) and macOS Sequoia with [`brew`](https://brew.sh). 

`MonopogenLite` was forked from `Monopogen` and edited as such to accommodate the work in [**MetaPlaq**](https://chanzuckerberg.com/science/programs-resources/cell-science/data-insights/metaplaq-integrative-single-cell-meta-analysis-for-atherosclerosis/). In **MetaPlaq** we meta-analyzed 140+ samples with single-cell RNA and ATAC sequencing data which have varying degrees of sequencing quality and depth. Many of these datasets are of limited use for reliable _de novo_ genotype calling. The main focus is on genetic ancestry inference and _cis_-acting expression quantitative trait loci (eQTL). 

`MonopogenLite` is a light-version of `Monopogen` and only includes the _germline calling_ and _phasing_ of single-nucleotide variants (SNVs) from single-cell sequencing data. 

It has a few improvements:

* ✅ Works cross-platform on macOS Sequoia and Linux Rocky8.
* ✅ Works with the newest versions of `bcftools`, `samtools`, `vcftools`, `tabix`, and `htslib`.
* ✅ Works with the version 4.1 of `BEAGLE` which still includes the `gl` option.
* ✅ Works with `python 3.9+`.
* ✅ Reproducible workflow to include a genome reference through [`refgenie`](http://refgenie.databio.org/en/latest/).
* ✅ Reproducible workflow to include 1000G phase 3 high-coverage b38 data including 3,202 individuals.
* ✅ Streamline code.
* ✅ Added `--debug`, `--logfile`, `--version` flags; messages from all steps, including parallel workers, go to one log file.
* ✅ Removed hard-coding of `--platform-library`; now works with _smartseq2_ and _celseq2_, aside of _10x_ (RNA and ATAC) data.
* ✅ Removed hard-coding of `--nthreads` for the phasing step executed through `BEAGLE`.
* ✅ One genotype per individual: every filtered BAM gets a single read group named after its sample ID, so samples are never merged and cells never become samples.
* ✅ Germline calling per chromosome with all samples of a study together; the per-chromosome BAM lists are built automatically from the `preProcess` output.
* ✅ Run germline calling step by step (`--step varScan`, `varProb`, `varPhasing`, or `all`); failed `bcftools` or `BEAGLE` steps are reported as errors.
* ✅ Handles STAR-based BAMs (`nM` tag), unmapped mates (e.g. 10x ATAC), and BAMs with `1` or `chr1` contig names.
* ✅ Optional UMI collapsing per cell (`--umi-collapse`, 10x RNA only; not recommended for germline calling).
* ✅ Added script to prepare outputs for TOPMed imputation.
* ✅ Code improved to only calls bi-allelic SNVs; no structural, INDELs or multi-allelic variants are called.
* ✅ Added handling of chromosome X - Phased haploid male samples in the phased imputation panel are now converted from "0" or "1" to "0|0" or "1|1", respectively.
* ✅ Refactored code to improve readability and maintainability, and cross-platform compatibility.
* ✅ Utilities work with regular `cpu` or `gpu`.
* ✅ Various utility-scrips:
  * ✅ `compare_vcf2ref.py` -- to compare a given VCF-file to a common reference VCF-file and list the non-overlapping variants.
  * ✅ `variant_ref_create.sh` -- create a reference VCF-file from a given set of BAM-files.
  * ✅ `variant_ref_checker.py` -- apply `bcftools stat` to calculate statistics on the newly generated reference file (including plotting). Including `variant_ref_checker_submit.sh` to submit the `variant_ref_checker.py` script to a job scheduler.
  * ✅ `makediploidmalesX.py` -- make haploid males diploid, or reverse male diploid genotypes to haploid on chromosome X in a given VCF-file. Including `makediploidmalesX_submit.sh` to submit the `makediploidmalesX.py` script to a job scheduler.
  * ✅ `cellsnp_vcf_merger.py` -- merge per-sample VCF files into 1 VCF per study. Including `cellsnp_vcf_merger_submit.sh` to submit the `cellsnp_vcf_merger.py` script to a job scheduler.

That said, `MonopogenLite` is just a light-weight version of `Monopogen`, originally developed by the [Ken Chen Lab](https://www.mdanderson.org/research/departments-labs-institutes/labs/ken-chen-laboratory.html) at the [MD Anderson Cancer Center](https://www.mdanderson.org/). All credit goes to the original authors and contributors of the concepts underlying `Monopogen`; we simply adapted and improved it in some functional and coding areas.

## 📦 Installation

`MonopogenLite` runs in one conda environment, `monopogen`, defined in [`environment.yml`](environment.yml). It holds the core tools (`python` 3.9, `pysam`, `samtools`, `bcftools`, `htslib`, `vcftools`, `openjdk`), the extras for the utility scripts (`matplotlib`, `cmcrameri`, `refgenie`) and `snakemake` 7.32 for the **MetaPlaq** workflows. `BEAGLE` 4.1 is included in the repository as `apps/beagle.jar`.

### 1. Get the code

```bash
git clone https://github.com/swvanderlaan/MonopogenLite.git
cd MonopogenLite
git checkout dev
```

### 2. Create the environment

With [`mamba`](https://github.com/conda-forge/miniforge) (or `conda`) on Linux or an Intel Mac:

```bash
mamba env create -f environment.yml
```

On a Mac with Apple Silicon (M1–M4) the same command creates a native (`osx-arm64`) environment; all packages are available for it. Only if you need an Intel (`osx-64`) environment, for example to match an existing one, create it this way; it then runs through Rosetta:

```bash
CONDA_SUBDIR=osx-64 mamba env create -f environment.yml
mamba activate monopogen
conda config --env --set subdir osx-64
```

The last line makes sure later installs into this environment are Intel packages too. Do not mix the two: packages built for Apple Silicon do not load in an Intel environment.

Tested package versions differ a little per platform: on Linux the solver picks `samtools`/`htslib` 1.20 and `pysam` 0.22.1 (the versions `MonopogenLite` was tested with), on macOS currently 1.25 and 0.24.1.

On an HPC cluster with `micromamba` (e.g. UVA Rivanna); the run scripts activate the environment by the name `monopogen`:

```bash
micromamba env create -f environment.yml
```

### 3. Check the installation

```bash
mamba activate monopogen
python src/MonopogenLite.py --version
samtools --version | head -1
bcftools --version | head -1
java -jar apps/beagle.jar 2>&1 | head -1
snakemake --version
```

### 4. Reference data

`germline` needs two resources:

- **Genome reference:** the GRCh38 FASTA the BAM files were aligned to (with its `.fai` index), e.g. `refdata-gex-GRCh38-2024-A/fasta/genome.fa` from 10x Genomics, or `refgenie seek hg38/fasta`.
- **Imputation panel:** a directory with one VCF per chromosome, named `1kGP_high_coverage_Illumina.SNVonly_poly.filtered_AF_5e-04.norm.fixvariantid.chrN.vcf.gz` (with `.tbi`) for chromosomes 1–22 and X. `src/variant_ref_create.sh` makes these as SLURM jobs: it downloads the phased high-coverage 1000 Genomes data (3,202 individuals), keeps bi-allelic SNVs with AF ≥ 0.0005 (`--af`, `--af-field`), normalises them against the reference (`.norm`), sets every variant ID to `CHROM:POS:REF:ALT` (`.fixvariantid`, the same IDs `germline` gives its calls), and also writes a combined sites-only panel for cellSNP (`...fixvariantid.chr1_22X.cellsnp.vcf.gz`). On Rivanna the finished panel is in `/project/cphg-millerlab/software/MonopogenLite/resources/`.

### Updating or repairing the environment

Install or change packages with `mamba`, not `pip`: `pip` can bring in versions that do not match (Snakemake 7 fails with `pulp` 2.8 or newer) or, on Apple Silicon, packages built for the wrong architecture. To bring an existing environment in line with `environment.yml`:

```bash
mamba env update -n monopogen -f environment.yml --prune
```

Python packages in your personal `~/.local` folder are visible in every environment and can override the environment's own. If something behaves differently for you than for others, check with `python -m pip list --user`, or run with `PYTHONNOUSERSITE=1`.

## 🚀 Usage

`MonopogenLite` runs in two steps: `preProcess` once per sample, then `germline` once per chromosome for all samples of a study together. See `python src/MonopogenLite.py preProcess --help` and `python src/MonopogenLite.py germline --help` for all options.

### 1. `preProcess` (per sample)

```bash
python src/MonopogenLite.py preProcess \
    --bamFile sampleA.bam.lst \
    --out study/monopogen_sampleA \
    --app-path apps \
    --platform-library 10x \
    --max-mismatch 3 \
    --nthreads 8
```

`--bamFile` is a comma-separated file with one line per sample: the sample ID and the path to a sorted, indexed BAM (`file.bam.bai` or `file.bai`). The sample ID becomes the sample name in the VCF.

```text
sampleA,/path/to/sampleA/outs/possorted_genome_bam.bam
```

For chromosomes 1–22 and X, `preProcess` writes `OUT/Bam/<sampleID>_chrN.filter.bam` (indexed). Reads are removed when they are unmapped, have `--max-mismatch` or more mismatches (`NM` tag, or `nM` for STAR-based BAMs; the default of 3 keeps reads with at most 2 mismatches, as in the original `Monopogen`), or are shorter than `--min-read-length` (default 30). Duplicate reads are kept: they add depth, and base and mapping quality filters in `bcftools` handle PCR errors.

| Data | Aligner | `--platform-library` | UMI | `--umi-collapse` |
| --- | --- | --- | --- | --- |
| 10x RNA-seq | CellRanger | `10x` | `UB` tag | allowed, not recommended |
| 10x ATAC-seq | CellRanger ATAC | `10x` | none | omit |
| SmartSeq2 | STARsolo | `smartseq2` | none | refused |
| CEL-Seq2 | e.g. bowtie2 | `celseq2` | in the read name | refused |

`--platform-library` does not change the filtering; it sets the read group fields of the filtered BAMs and whether `--umi-collapse` is allowed.

### 2. `germline` (per chromosome, all samples of a study)

```bash
python src/MonopogenLite.py germline \
    --chromosome 1 \
    --region region.lst \
    --out study \
    --reference genome.fa \
    --imputation-panel resources/ \
    --app-path apps \
    --nthreads 8
```

`--out` is the study folder holding the `monopogen_<sampleID>` folders from `preProcess` (or set `--preprocess-dir`). If `OUT/ListedBamsPerChr/chrN.lst` does not exist yet, `germline` builds it from all `*/Bam/*_chrN.filter.bam` files it finds there.

The region file has one region per line: a chromosome (whole chromosome) or `chromosome,start,end`. Chromosomes may be written with or without `chr` and in any order; a chromosome may have several regions, and lines starting with `#` are ignored.

```text
chr1
chr2,1,120000000
chr2,120000001,242193529
X
```

For each region, `germline` calls bi-allelic SNVs with `bcftools mpileup` and `call` (mapping and base quality ≥ 20), then estimates genotype probabilities and phases with `BEAGLE` against the imputation panel. `--step varScan`, `varProb` or `varPhasing` runs a single step; the default is `all`. Output in `OUT/`:

| File | Content |
| --- | --- |
| `monopogen_<region>_sorted.vcf.gz` (+ `.tbi`) | Unphased calls, all samples |
| `monopogen_<region>.gp.vcf.gz` | Genotype probabilities (`BEAGLE`) |
| `monopogen_<region>.phased.vcf.gz` | Phased genotypes (`BEAGLE`) |
| `runGermline_<region>.sh` | The job script that produced them |

A region without any called variant skips `BEAGLE` and keeps only the unphased VCF.

### Running on a SLURM cluster

`src/runPreprocess.sh` runs `preProcess` as a SLURM array with one task per sample; `src/runGermline.sh` runs `germline` as an array with one task per chromosome (`--array=1-23`, where 23 is X). Both use `PROJECT_DIR/<study>/monopogen_<sample>` as the folder layout. Adjust the `#SBATCH` lines and default paths for your cluster.

```bash
sbatch src/runPreprocess.sh --input-dir /path/to/cellranger_outs/<study> --study <study> --platform 10x
sbatch src/runGermline.sh --study <study>
```

## 📣 Acknowledgements
We thank Damian Ulrich for the per-chromosome germline calling and phasing, and Teresa Jhee for the per-chromosome read counting and contig-name handling in `preProcess`; both also contributed to testing `MonopogenLite` on the **MetaPlaq** studies.

Dr. Sander W. van der Laan is funded through EU H2020 TO_AITION (grant number: 848146), EU HORIZON NextGen (grant number: 101136962), EU HORIZON MIRACLE (grant number: 101115381), Health~Holland PPP Allowance ‘Getting the Perfect Image’, and CZI ['MetaPlaq'](https://chanzuckerberg.com/science/programs-resources/cell-science/data-insights/metaplaq-integrative-single-cell-meta-analysis-for-atherosclerosis/).

We are thankful for the support of the Leducq Fondation ‘PlaqOmics’ and 'AtheroGen'. The research for this contribution was made possible by the AI for Health working group of the [EWUU alliance](https://aiforhealth.ewuu.nl/). The collaborative project ‘Getting the Perfect Image’ was co-financed through use of PPP Allowance awarded by Health~Holland, Top Sector Life Sciences & Health, to stimulate public-private partnerships. Part of the work and data generation (scRNAseq and bulk RNAseq) were funded through ERA-CVD 'druggable-MI-targets' project (grantnumber: 01KL1802). 

Plaque samples are derived from carotid endarterectomies as part of the [Athero-Express Biobank Study](https://www.atheroexpress.nl)([Verhoeven _et al._](https://doi.org/10.1007/s10564-004-2304-6)) which is an ongoing study in the UMC Utrecht.

### 📝 Disclosures
Dr. Sander W. van der Laan has received Roche funding for unrelated work.

<a href='https://uefconnect.uef.fi/en/group/miracle/'><img src='images/UEF_Miracle_Logo-07.png' align="center" height="75" /></a> <a href='https://www.to-aition.eu'><img src='images/to_aition.png' align="center" height="75" /></a> <a href='https://www.health-holland.com'><img src='images/logo_NL_HealthHollland_Wit-Oranje_RGB.png' align="center" height="35" /></a> <a href='https://www.nextgentools.eu'><img src='images/NextGen_1_Red.png' align="center" height="35" /></a> <a href='https://www.era-cvd.eu'><img src='images/ERA_CVD_Logo_CMYK.png' align="center" height="75" /></a> <a href=''><img src='images/leducq-logo-large.png' align="center" height="75" /></a> <a href='https://www.fondationleducq.org'><img src='images/leducq-logo-small.png' align="center" height="75" /></a> <a href='https://osf.io/zcvbs/'><img src='images/worcs_icon.png' align="center" height="75" /></a> <a href='https://doi.org/10.1007/s10564-004-2304-6'><img src='images/AE_Genomics_2010.png' align="center" height="100" /></a>

#### 📦 Changes log
    
    _Version:_      v1.4.11</br>
    _Last update:_  2026-11-09</br>
    _Changes:_      see CHANGELOG.md</br>
    _Written by:_   Jinzhuang Dou | jdou1 [at] mdanderson [dot] org; Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.</br>
    _Contributors:_ Damian Ulrich; Teresa Jhee.
    
    **MoSCoW To-Do List**
    The things we Must, Should, Could, and Would have given the time we have.
    _M_
    - [ ] test v1.4.9 on Rivanna with a 10x ATAC study (Turner) and a SmartSeq2 study (Katyayani)

    _S_
    - [x] remove the hidden `--max-softClipped` option once no script passes it anymore (v1.4.11)
    
    _C_

    _W_
    - [ ] add the support for somatic SNV calling and LD refinement  

--------------

#### MIT License (MIT)
##### Copyright (c) 1979-2026 Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.
