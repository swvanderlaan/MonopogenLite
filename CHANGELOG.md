# 📜 CHANGELOG.md — MonopogenLite

All notable changes to `MonopogenLite` will be documented in this file.


## 🛠️ v1.3.3 -- 2025-06-10
### Smaller Fixes
- 🐛 `MonopogenLite.py`: Fixed `NameError` in `germline()` by adding missing `samples = read_sample_list_file(args.region)` before calling `build_sample_commands()`.
- ✅ Patch added: defined `read_sample_list_file()` to resolve the `NameError`. This function extracts sample IDs from the `--region` file. 


## 🛠️ v1.3.2 -- 2025-06-07
### Smaller Fixes
- 🐛 `MonopogenLite.py`: Fixed `FileNotFoundError` by ensuring log file directory is created before initializing the logger.


## 🛠️ v1.3.1 -- 2025-06-05
### Improved doc strings
- 📝 Improved doc strings in `MonopogenLite.py` for better clarity and understanding of functions.
### Smaller Fixes
- 🐛 Fixed issue with indentations. 


## 🛠️ v1.3.0 -- Major Patch -— 2025-06-05
### Major Enhancements & Refactors in `MonopogenLite.py`
- 🧱 Refactored command generation into helper functions:
- `generate_bcftools_command()`, `generate_beagle_cmd_gp()`, and `generate_beagle_cmd_gt()`.
- 📁 Created utility `prepare_output_dirs()` for consistent directory creation in both `germline()` and `preProcess()`.
- 📁 `preProcess()` now uses shared `prepare_output_dirs()` function for cleaner and more consistent output directory setup.
- 📉 Removed global dependency on `out`; now passed explicitly to improve modularity.
- 🧪 Refactored `germline()` and `preProcess()` logic by introducing two new helper functions: 
  - `build_sample_commands()` for modular command generation
  - `write_job_script()` for dynamic SLURM script creation. 
  This significantly improves clarity, maintainability, and unit test potential.
- 📝 Refactored SLURM script generation to use Python `f-strings` for dynamic version and job ID injection (e.g., `v{VERSION}` and `{jobid}`).
- 🧰 Introduced helper function `write_slurm_script()` to centralize and standardize SLURM job script creation across `germline()` and `preProcess()` workflows.
### Improved
- 🧠 **`germline.py`**: Enhanced validation error messages and traceback logging.
- 🧪 **`germline.py`**: Refactored read group header assignment: replaced float `LB=0.1` with string-based sampleID and added fallback for unknown platforms.
#### Improved Logging
- 🧾 Added `--logfile` argument to write logs to a file (default: `OUTDIR/MonopogenLite.YYYYMMDD.log`).
- 🧪 Appended log file output to console when `--verbose` is used.
- 🔧 Refactored `setup_logger()` to configure both console and file logging consistently.
- 🔊 Moved --logfile handling to the top-level `main()` function so that both germline() and preProcess() write to the same log file consistently.
- 🧪 Logged versions of key tools (`bcftools`, `beagle.jar`) at runtime when `--verbose` is enabled, enhancing reproducibility.
#### Better Error Handling
- 🛑 Replaced fragile `assert` statements with explicit validation and informative `sys.exit(1)` errors.
- 📦 Validated presence of `beagle.jar` and provided a clear error message if not found.
#### Code Quality Fixes
- 🧾 Standardized all path constructions using `os.path.join()` instead of `+ '/' +` for better reliability and cross-platform compatibility.
- 🧱 Refactored all file and directory operations to use `pathlib.Path` instead of `os.path` and string concatenation, improving cross-platform reliability and code readability. Standardized path variable naming to out_path for consistency.
- ✅ Removed unnecessary parentheses in if conditions (`if(len(record))` → if `len(record))`.
- 🧮 Replaced inefficient list comprehension with a set() for chromosome checks (in `[f"chr{n}"...]` → in `{f"chr{n}"...})`.
#### 📜 Shell Script Improvements
- 📝 All generated job scripts now include:

    ```
    #!/bin/bash
    echo "Running: $(basename $0)"
    echo "Start time: $(date)"
    ```
### Smaller Fixes
- 🐛 Fixed issue in `runGermline.sh` where double quotes would throw errors or silently misinterpret the python command.
- 🐛 **`germline.py`**: Fixed typo in `robust_get_tag()` where `tagname` was undefined.
- 🔍 **`germline.py`**: Improved logic in `BamFilter()` to safely strip `'chr'` prefix only when needed.
- 🔧 **`germline.py`**: Updated `validate_sample_list_file()` to raise `ValueError` instead of `sys.exit(1)` for better error propagation.
- 📛 **`germline.py`**: Prevented double `'chr'` prefix in `addChr()` header rewriting.


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