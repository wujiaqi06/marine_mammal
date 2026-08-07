# Marine mammal genome evolution: NC reproducibility package

This repository contains the custom source code, lightweight source data,
documentation, branch-label metadata and validation checks for the marine
mammal genome-evolution manuscript submitted to *Nature Communications*.
Release `v1.1-nc-submission` extends the frozen `v1.0-endpoint-fix` analysis
baseline with the current NC Fig. 6 species/ancestral projections,
Supplementary Fig. S3 ancestral fingerprints and reviewer-facing
reproducibility checks. Scientific values from the endpoint-fix analysis are unchanged.

Large study-generated inputs and outputs are deposited in Dryad under DOI
<https://doi.org/10.5061/dryad.dz08kpsd4>. The DOI is reserved and public
release is pending manuscript publication. During peer review, reviewers
should use the private Dryad link provided in the manuscript. The private
reviewer URL is intentionally not stored in this public repository.

SplitAligner is maintained separately at
<https://github.com/wujiaqi06/SplitAligner>.

## Reviewer quick start

The bundled demo uses real lightweight Source Data and does not require the
large Dryad matrices.

```bash
# From GitHub:
git clone --branch v1.1-nc-submission --depth 1 \
  https://github.com/wujiaqi06/marine_mammal.git
cd marine_mammal
bash demo/run_demo.sh

# Or, from the journal editor/reviewer ZIP, enter its marine_mammal directory:
cd marine_mammal
bash demo/run_demo.sh
```

Expected final lines:

```text
PASS NC source-data checks
PASS demo output matches locked expected output
```

The demo verifies the main endpoint-fix source tables and the current NC
Fig. 6/Fig. S3 source tables. On the tested desktop it completes in less than
5 seconds after R is available. See `demo/README.md` for the generated output.

The public repository is <https://github.com/wujiaqi06/marine_mammal>.
Tag `v1.0-endpoint-fix` preserves the endpoint-fix baseline. Tag
`v1.1-nc-submission` adds the current Fig. 6 and Supplementary Fig. S3 code,
lightweight Source Data and reviewer-facing reproducibility checks without
changing the endpoint-fix scientific outputs.

## System requirements

### Lightweight demo and source-data checks

- macOS or Linux with a POSIX shell.
- R 4.2 or later; tested with R 4.4.2 on macOS Sonoma 14.4 (Apple silicon).
- No non-standard hardware. One CPU core and less than 1 GB RAM are sufficient.

### Figure regeneration and downstream analyses

- R 4.4.2, with tested packages listed in `env/R_required_packages.tsv`.
- Python 3.9.6, with tested packages listed in
  `env/python_requirements.txt`.
- Perl 5.34.1 for alignment-QC and SplitAligner scripts; the main SplitAligner
  workflow uses core Perl modules only.
- Upstream phylogenomic regeneration additionally requires IQ-TREE2 2.1.3,
  ASTRAL III and PAML/BaseML 4.9j as described in Methods and
  `docs/reproduction_guide.md`.

The full upstream tree, BaseML, SplitAligner and nested-LASSO workflows are
compute- and storage-intensive. They are not the desktop demo: use a multicore
Unix-like workstation or cluster and allow workflow-specific hours to days.
The submitted analyses used up to 12 CPU cores where parallel execution was
implemented.

## Installation

The code is script based; no compilation is required for the bundled demo.

```bash
# Install missing R dependencies used by downstream analysis/figure scripts.
Rscript env/install_R_dependencies.R

# Optional isolated Python environment for figure scripts.
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r env/python_requirements.txt
```

Repository checkout normally takes less than 1 minute. Installing the R and
Python dependencies typically takes approximately 5-15 minutes on a normal
desktop with a working package mirror and network connection; this varies by
platform and whether binary packages are available.

Exact environment captures from package preparation are retained under
`env/`. They provide provenance and are not a claim that every upstream
phylogenomic tool was run in that same local session.

If Matplotlib cannot write its normal font cache on a restricted system, use a
writable cache directory, for example:

```bash
mkdir -p /tmp/marine_mammal_mpl
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/marine_mammal_mpl \
  python3 code/12_NC_species_ancestral_fingerprints/build_Figure6_layout.py
```

## Demonstration dataset and expected output

`source_data/` contains small real tables supporting the manuscript figures.
The NC additions are:

- `source_data/Fig6_species_ancestral_fingerprints/`
- `source_data/FigS3_ancestor_fingerprints/`

`bash demo/run_demo.sh` writes `demo/output/demo_summary.tsv` and compares it
with `demo/expected_output/demo_summary.tsv`. Expected counts include six
terminal species decomposed in Fig. 6B/C, 71 marine predictors, 101 binary
aquatic-dependence predictors and 13 ancestral profiles in Fig. S3.

## Reproducing manuscript analyses

1. Download and unpack the paired Dryad dataset.
2. Set the endpoint-fix project root used by full-workflow scripts:

   ```bash
   export MARINE_MAMMAL_ENDPOINTFIX_ROOT=/path/to/unpacked_endpointfix_project
   export MARINE_MAMMAL_NC_OUTPUT_ROOT=/path/to/reproduction_outputs
   ```

3. Follow `docs/reproduction_guide.md` and `docs/reproduce_map.tsv` for the
   stage-specific command and source/output mapping.
4. Run the package checks:

   ```bash
   bash tests/smoke_tests/run_all.sh
   ```

The full scripts expect the endpoint-fix directory structure documented in
the reproduction guide. `config/example_paths.yaml` is a human-readable path
template only; scripts use the environment variables documented above rather
than parsing that YAML file.

## Using the code with other data

SplitAligner accepts a fixed species tree and gene trees as documented in its
own repository and in `code/04_SplitAligner/README.md`. Downstream scripts
require gene-by-branch matrices whose branch identifiers follow the endpoint-
fix mapping contract, trait tables with the documented columns, and terminal
species labels that match the tree. Users adapting the downstream scripts
should preserve these schemas and update the environment-controlled roots.

## Branch-label mapping warning

SplitAligner's branch-label annotation rule changed after the downstream
manuscript analyses were built. Current direct SplitAligner outputs must be
mapped to the old downstream branch-label rule/order before reproducing the
t-test, deterministic ASR joins, LASSO, Source Data or figures. Dryad branch-
coordinate and GBI matrices are already converted to the old labels/order.
Crosswalks, a portable converter and a 601-branch smoke test are provided in
`metadata/branch_label_mapping/`. Example:

```bash
Rscript metadata/branch_label_mapping/apply_branch_label_crosswalk.R \
  --input=/path/to/current_split_aligner_matrix.tsv \
  --output=/path/to/old_label_matrix.tsv
```

Do not convert the paired Dryad matrices again.

## Evidence boundaries for Fig. 6 and Fig. S3

Fig. 6 terminal and internal-branch coordinates are descriptive fitted
projections onto the final full-data fingerprints after validation. Fig. 6B/C
and Fig. S3 contributions are scaled GBI multiplied by the corresponding
full-data coefficient; the intercept is omitted from bars but included in
profile-score reconstruction. Internal-branch profiles are not direct ancestral
habitat assignments. Final panel typography was adjusted in Adobe Illustrator;
the included scripts regenerate source-consistent draft artwork and all plotted
scientific values.

## Key endpoint-fix values

- 302 terminal species.
- Fig. 2 nested AUCs: marine 0.936; binary aquatic-dependence 0.826.
- Final full-data fingerprints: 71 marine and 101 binary aquatic-dependence
  predictors, with 148 genes in their union.
- Fig. 5B positive-count-matched permutation: 894 FDR genes, 98.0% slow.

## License and external-data boundary

Custom code written for this study is released under the MIT License; see
`LICENSE`. External TOGA-derived mammalian coding-gene alignments from Michael
Hiller Lab resources are not redistributed here and are not relicensed by this
repository. Study-generated large processed data products are provided through
the paired Dryad dataset.

## Citation

Please cite the associated manuscript, the Dryad dataset
<https://doi.org/10.5061/dryad.dz08kpsd4>, this code release and SplitAligner.
Machine-readable citation metadata are in `CITATION.cff`.
