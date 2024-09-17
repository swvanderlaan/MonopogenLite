# MonopogenLite
Germline SNV calling and phasing from single-cell sequencing data (for macOS Sequoia and linux Rocky8).


This is a fork of the original [`MonopogenLite`](https://github.com/KChen-lab/Monopogen) which works with python (3.7+), and `samtools`, `vcftools`, and `bcftools`, as well as in the context of `Rocky8` (Linux) and macOS Sequoia with [`brew`](https://brew.sh). 

`MonopogenLite` was forked from `Monopogen` and edited as such to accommodate the work in [**MetaPlaq**](https://chanzuckerberg.com/science/programs-resources/cell-science/data-insights/metaplaq-integrative-single-cell-meta-analysis-for-atherosclerosis/). In **MetaPlaq** we meta-analyzed 140+ samples with single-cell RNA and ATAC sequencing data which have varying degrees of sequencing quality and depth. The main focus is on genetic ancestry inference and _cis_-acting expression quantitative trait loci (eQTL). Thus, `MonopogenLite` is a light-version of `Monopogen` and only includes the germline calling and phasing of single-nucleotide variants (SNVs) from single-cell sequencing data. It has a few improvements:

* [x] works with the newest versions of `bcftools`, `samtools`, `vcftools`, and `htslib`
* [x] works with the version 4.1 of `beagle` which still includes the `gl=` option
* [x] works with `python 3.9+`
* [x] reproducible workflow to include a genome reference through [`refgenie`](http://refgenie.databio.org/en/latest/) 
* [x] reproducible workflow to include 1000G phase 3 high-coverage b38 data including 3,202 individuals
* [x] streamline code 
* [x] added `--debug` and `--version` flags
* [x] removed hard-coding of `--platform-library`; now works with _smartseq2_ and _celseq2_, aside of _10x_ data
* [x] removed hard-coding of `--nthreads` for the phasing step executed through `beagle`
* [x] improved handling of per-sample or per-cell SNV calling and phasing
* [x] added script to prepare outputs for TOPMed imputation
* [x] code improved to only calls bi-allelic SNVs, no structural, INDELs or multi-allelic variants are called
* [ ] added handling of chromosome X - Issue to tackle: "To ensure proper ploidy of male samples in the phased panel, we converted “0|1”, “1|0”, and “1|1” GTs into a haploid representation (i.e. “1”) in nonPAR regions of chrX in males in the new (v2) version of the chrX VCF." We need to somehow filter out those haploid representations from the nonPAR region - in other words ignore the nonPAR region (for now). 

# Acknowledgements
Dr. Sander W. van der Laan is funded through EU H2020 TO_AITION (grant number: 848146), EU HORIZON NextGen (grant number: 101136962), EU HORIZON MIRACLE (grant number: 101115381), Health~Holland PPP Allowance ‘Getting the Perfect Image’, and CZI ['MetaPlaq'](https://chanzuckerberg.com/science/programs-resources/cell-science/data-insights/metaplaq-integrative-single-cell-meta-analysis-for-atherosclerosis/).

We are thankful for the support of the Leducq Fondation ‘PlaqOmics’. The research for this contribution was made possible by the AI for Health working group of the [EWUU alliance](https://aiforhealth.ewuu.nl/). The collaborative project ‘Getting the Perfect Image’ was co-financed through use of PPP Allowance awarded by Health~Holland, Top Sector Life Sciences & Health, to stimulate public-private partnerships. Part of the work and data generation (scRNAseq and bulk RNAseq) were funded through ERA-CVD 'druggable-MI-targets' project (grantnumber: 01KL1802). 

Plaque samples are derived from carotid endarterectomies as part of the [Athero-Express Biobank Study](https://doi.org/10.1007/s10564-004-2304-6) which is an ongoing study in the UMC Utrecht.

## Disclosures
Dr. Sander W. van der Laan has received Roche funding for unrelated work.

<a href='https://www.era-cvd.eu'><img src='images/ERA_CVD_Logo_CMYK.png' align="center" height="75" /></a> <a href='https://www.plaqomics.com'><img src='images/leducq-logo-large.png' align="center" height="75" /></a> <a href='https://www.fondationleducq.org'><img src='images/leducq-logo-small.png' align="center" height="75" /></a> <a href='https://osf.io/zcvbs/'><img src='images/worcs_icon.png' align="center" height="75" /></a> <a href='https://doi.org/10.1007/s10564-004-2304-6'><img src='images/AE_Genomics_2010.png' align="center" height="100" /></a>

#### Changes log
    
    _Version:_      v1.0.0</br>
    _Last update:_  2024-09-17</br>
    _Written by:_   Jinzhuang Dou | jdou1 [at] mdanderson [dot] org; Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.
    
    **MoSCoW To-Do List**
    The things we Must, Should, Could, and Would have given the time we have.
    _M_
    
    - [] add support for chromosome X
    
    _S_

    - [] add the support for somatic SNV calling and LD refinement
    
    _C_

    _W_

    **Changes log**
    * v1.0.0 Initial version. 

--------------

#### Creative Commons BY-NC-ND 4.0
##### Copyright (c) 1979-2024 Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.
