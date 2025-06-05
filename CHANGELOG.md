# 📜 CHANGELOG.md — MonopogenLite

All notable changes to MonopogenLite will be documented in this file.

## 🛠️ v1.2.6 — 2025-06-05
### Fixed
- 🐛 Fixed `runPreprocess.sh` and `runGermline.sh` to be more flexible, including submission-scripts.
### Updated
- 🏷️ Updated file names of submission scripts.
### Added
- 🧬 Added `cellsnp_vcf_merger.py` (and associated `cellsnp_vcf_merger_submit.sh`) to merge per-sample VCF files into 1 VCF per study.

## 🛠️ v1.2.5 — 2025-06-03
### Fixed
- 🐛 **`germline.py`**: Fixed an issue where the variable `val` could be referenced before assignment in the `BamFilter` function. Now properly checks `val is not None` before comparing mismatch values.
### Added
- 🏷️ Added `CHANGELOG.md` for better tracking of changes.
### Updated
- 🏷️ Updated file names of submission scripts.

## 🧭 v1.2.4 — 2024-10-01
### Fixed
- 🔧 **`MonopogenLite.py`**: Fixed the hardcoded reference to the imputation panel.

## 📘 v1.2.3 — 2024-09-19
### Added
- 📣 **`runGermline.sh`**: Added more informative logging and output to the `runGermline` scripts.

## 🧬 v1.2.2 — 2024-09-18
### Updated
- 🏷️ Updated `germline.py` to account for `RG` header tags based on the sequencing platform (e.g. 10x, SmartSeq2, CEL-Seq2).

## 🧮 v1.2.1 — 2024-09-18
### Added
- 📊 Added a script (`compare_vcf2ref.py`) to count overlapping variants between BAM-derived inputs and resulting VCF outputs.

## 🚀 v1.2.0 — 2024-09-18
### Added **`MonopogenLite.py`**: 
- 📏 Support for minimum read length (default: 30 bp).
- 🔬 Support for specifying sequencing platform (`10x`, `smartseq2`, `celseq2`).
- 🧪 UMI collapsing support via UMI tags.

## 🧬 v1.1.0 — 2024-09-18
### Added **`MonopogenLite.py`**: 
- 🧬 Chromosome X support added for variant processing. These are `makediploidmalesX.py` and `makediploidmalesX_submit.sh`.
- 🧪 Added script to create a reference (`variant_ref_create.sh`).
- 🧪 Added script to check the reference created (`variant_ref_checker.py` and `variant_ref_checker_submit.sh`).

## 🎉 v1.0.0 — 2024-09-17
### Initial Release
- 🌱 **`MonopogenLite.py`**: Initial release of **MonopogenLite**, a lightweight fork of Monopogen focused exclusively on germline variant calling.