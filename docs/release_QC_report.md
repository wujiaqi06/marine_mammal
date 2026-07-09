# Release QC report for v1.0-endpoint-fix

Generated: 2026-07-09 11:05:52 JST

## Summary

- Private Dryad reviewer URL: not present.
- Local absolute `/Users/...` paths: not present outside this QC report.
- Files larger than 100 MB outside `.git`: none.
- Core nested gLOOCV, corrected LASSO, fold-wise ASR sensitivity and figure-generation scripts are included under `code/` with path roots controlled by environment variables or local config.
- Non-portable provenance scripts with local paths were excluded from the GitHub release; `provenance_scripts/nonportable_script_manifest.tsv` records the excluded provenance-only scripts.
- Smoke tests: source-data presence PASS; endpoint-fix value check PASS; forbidden legacy-value check PASS.
- Remaining grep hits are false positives from biological terms (`secretion`, `secreted`), code variable names (`token`) or `.gitignore` safety patterns.

## Required private/safety grep

Command:
```bash
grep -RInE "LINK_NOT_FOR_PUBLICATION|private reviewer|reviewer link|cover letter|suggested reviewers|submitted elsewhere|desk reject|rejection|Science submission|Nature submission|NEE|Claude|ChatGPT|OpenAI|Anthropic|api[_-]?key|token|password|secret" . --exclude-dir=.git --exclude="release_QC_report.md" || true
```

Output / interpretation:
```text
./source_data/Fig3_marine_single_gene_screen/endpointfix_Fig3C_enrichment_terms_for_plot.tsv:57:marine_slow	KEGG	hsa04927	Cortisol synthesis and secretion	15	65	0.52	0.36	0.0274	1.562249437179612	False	excluded_from_main_figure	Select STRING/Monarch phenotype and measurement terms, allowing EFO and HP ontology terms; prioritize hematological, erythrocyte, hemoglobin, platelet, blood-cell, cardiovascular, anthropometric/body-measure, or protein-measurement wording; exclude broad/generic terms.
./source_data/Fig4_corrected_full_data_architecture/Figure4C_TableS5_annotation_table_final.tsv:59:ASMTL	shared	True	True	-0.225475	-0.00682425	0.225475	marine slow direction negative; aquatic slow direction negative	Unassigned / low-confidence	Endocrine / circadian / systemic regulation	ASMT/ASMTL-locus; melatonin-pathway-adjacent methyltransferase-like candidate	low	False	False	True	False	ASMTL is retained as an ASMT/ASMTL-locus and melatonin-pathway-adjacent context-only candidate. ASMTL ranks 262nd by marine t-test p-value, whereas ASMT ranks 3513th and is not significant in the current screen; therefore ASMTL should not be interpreted as direct evidence for altered melatonin secretion. Use only as cautious context in Supplementary Table S5 or Discussion.	NCBI Gene; HGNC; PMID:39859321; NCBI Gene; HGNC	ASMTL is retained as an ASMT/ASMTL-locus and melatonin-pathway-adjacent context-only candidate. ASMTL ranks 262nd by marine t-test p-value, whereas ASMT ranks 3513th and is not significant in the current screen; therefore ASMTL should not be interpreted as direct evidence for altered melatonin secretion. Use only as cautious context in Supplementary Table S5 or Discussion.	ASMTL is retained as an ASMT/ASMTL-locus and melatonin-pathway-adjacent context-only candidate. ASMTL ranks 262nd by marine t-test p-value, whereas ASMT ranks 3513th and is not significant in the current screen; therefore ASMTL should not be interpreted as direct evidence for altered melatonin secretion. Use only as cautious context in Supplementary Table S5 or Discussion.	Table S5-only ASMT/ASMTL-locus context; not displayed in main Figure4C and not evidence for direct altered melatonin secretion.
./source_data/Fig4_corrected_full_data_architecture/Figure4C_TableS5_annotation_table_final.tsv:62:SPINT4	aquatic_specific	False	True	0.0	-0.0769933	0.0769933	aquatic slow direction negative	Epithelial / body-surface interface	Epithelial / body-surface interface	secreted serine protease inhibitor / epithelial protease regulation	medium	True	True	False	False	secreted serine protease inhibitor / epithelial protease regulation. Serine protease inhibitor family member; plausible epithelial/protease-interface placement.	UniProt; NCBI Gene; GO	Serine protease inhibitor family member; plausible epithelial/protease-interface placement.	Manual review: limited functional literature; keep medium confidence.	Main Figure4C display candidate under compressed-module policy.
./source_data/Fig4_corrected_full_data_architecture/Figure4C_TableS5_annotation_table_final.tsv:90:APH1B	marine_specific	True	False	-0.152261	0.0	0.152261	marine slow direction negative	Membrane trafficking / ER / vesicle processing	Membrane trafficking / ER / vesicle processing	gamma-secretase complex subunit / membrane proteolysis	high	True	False	True	False	gamma-secretase complex subunit / membrane proteolysis. Membrane protein processing/gamma-secretase complex.	UniProt; NCBI Gene; Reactome	Membrane protein processing/gamma-secretase complex.		Optional/Table S5-only by default because membrane trafficking / ER / vesicle processing is excluded from main Figure4C.
./data_manifest/source_data_manifest.tsv:2:source_data/Baseline_slow_gene_layer_overlap/Supplementary_Table_S4_baseline_slow_gene_layer_overlap.tsv	5686	a69a22d57e3095c5e58f98744ff6bd130c86900985f7bfcf88ea8774f274faf7	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:3:source_data/Baseline_slow_gene_layer_overlap/slow_gene_direction_crosstab.tsv	7356	ab2c06b60dd08d1c10485e3fcaf05e8a1556213c3cff4562760166a73e8ba91b	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:4:source_data/Baseline_slow_gene_layer_overlap/slow_gene_overlap_partitions.tsv	134698	a35bb25cce49e4ea1a40ea730c6b697b77e083151722069352a81d4bcd6ad2e8	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:5:source_data/Baseline_slow_gene_layer_overlap/slow_gene_set_counts.tsv	1952	2a4b3ce9f0e58fa1dd1dec3ad0ad3f512cbef9eed9fb049c07514f17984ea68d	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:6:source_data/Fig1_trait_framework/Figure1C_score_clade_counts.tsv	434	3894729dacadc14dd28c9eef12fdc009273816bfa3ddba7d429aecf9eb0d0137	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:7:source_data/Fig1_trait_framework/trait_table.mammal302.active_TY_NK_final_18pt.tsv	20378	b7b04906d074bb2bd8f82b5b646fcd76b1b880e566ef6cb8a0cf0ef54ed0ccd2	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:8:source_data/Fig2_nested_gLOOCV/Figure2B_nested_ttest_ROC_data.tsv	177197	f99870f44aa2e0cd6863135de08f368194b317a702af6aeb1fffb38e2256b0b7	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:9:source_data/Fig2_nested_gLOOCV/Figure2C_nested_ttest_OOF_distribution_data.tsv	222948	4f4d4d59c311bbfa4802475eec6c1abdfb421cee0af7e7e1dd8075e2137c37bd	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:10:source_data/Fig2_nested_gLOOCV/Figure2_nested_ttest_plot_summary.tsv	1308	5bce1339897a8b1baf9f3d2ea31a1ab30448851181908b3bff4b6d4f818d99b1	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:11:source_data/Fig3_marine_single_gene_screen/endpointfix_Fig3A_screening_plot_table.tsv	4605635	68f30f12d1ea14bedb0b0e52bcc45e27cbf3be0324df60312ddde0826357386f	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:12:source_data/Fig3_marine_single_gene_screen/endpointfix_Fig3B_directional_asymmetry_summary.tsv	2563	0765a82dd8d7e83d8684d2ba87d034e0a6bc72cdc83f94902ed17d79144b2ee0	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:13:source_data/Fig3_marine_single_gene_screen/endpointfix_Fig3C_enrichment_terms_for_plot.tsv	61552	9c7f5198c45d9f40067f5ae27d38d78170a540a37a909da26af0d0b84283ffb7	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:14:source_data/Fig4_corrected_full_data_architecture/Figure4A_plot_table.tsv	215	d677cad7bb955e72f71c23fec989f267a0cd8f4bf0f4acdbaf4f831794f6ba2a	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:15:source_data/Fig4_corrected_full_data_architecture/Figure4B_plot_table.tsv	35425	5acf9a8772a2e30f4abf8c4b517aade8b1297431f2295c1af674a28fd55e4157	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:16:source_data/Fig4_corrected_full_data_architecture/Figure4C_TableS5_annotation_table_final.tsv	80849	217579c1d406e835919e44152f4516297d269a0c4842ac16b9ea066cce8cc8f0	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:17:source_data/Fig4_corrected_full_data_architecture/Figure4C_main_display_cell_summary_final.tsv	4459	ec4fead3b9a7ff11fa6a2e6237c9026c29c57efabca15ddec1dd54e79792363e	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:18:source_data/Fig4_corrected_full_data_architecture/Figure4C_plot_table.tsv	13619	4ae7735b2df9190bb757d28ab6f4b5a8fca4b48446c99f3199cc5c9ee05d35be	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:19:source_data/Fig5_sensitivity_permutation_turnover/Figure5A_nested_sensitivity_plotdata.tsv	3108	231383e8934765840b8b3780292597d39e0f641397e26bdab0c5ee678df1911c	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:20:source_data/Fig5_sensitivity_permutation_turnover/Figure5A_nested_sensitivity_table.tsv	2411	0dcd273bdfe26909430950d1aa66a917fecd40fdca31b1721bf21601a8450cad	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:21:source_data/Fig5_sensitivity_permutation_turnover/Figure5B_endpointfix_positive_count_permutation_plotdata.tsv	106064	bb00719279a55ea09e326feea8fb79cc6f3feb9c3b2b03ec3bf7e79b54f7a692	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:22:source_data/Fig5_sensitivity_permutation_turnover/Figure5B_endpointfix_source_check.tsv	1271	6d4307377b4f5ae083bb521cb8782d9631f292a2354f408950a76aa366b23e90	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:23:source_data/Fig5_sensitivity_permutation_turnover/Figure5C_cross_layer_comparison_supplementary.tsv	805	3b4332de8c74cdacddeb54b3d08431fb4704b84bcea088f55b06c6b81712f633	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:24:source_data/Fig5_sensitivity_permutation_turnover/Figure5C_module_null_summary_2col.tsv	1156	91b59411a85df88ebd66ef9885a974d8c614372b99f57458a23d39707e3f539e	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:25:source_data/Fig5_sensitivity_permutation_turnover/Figure5C_predictor_turnover_metrics_2col.tsv	478	cb4836c107569272bdd9ae8b5756e74ba500fe61ffb218afe82560be8f0f343a	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:26:source_data/FigS1_aquaticity_PCA/FigureS1_PCA_plot_tables.tsv	34602	5329e95740ef0f325a3982f5a6b197340ea3adb0b00cd0bc5058bd0ba452e082	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:27:source_data/FigS1_aquaticity_PCA/FigureS1_PCA_summary.txt	755	42a0a6c54a7329a6f2291c6fc4d6c1fa6ea8bc4bd1df22f19dce9947bfc558c7	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:28:source_data/FigS1_aquaticity_PCA/aquaticity_5d_scores.mammal302.slim_TY_NK_final_18pt.tsv	137354	bf96dcb73b6d66c89db1f25e5f77275852e627a4f133568fd5a19719630d0f47	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:29:source_data/FigS2_aquatic_single_gene_screen/endpointfix_FigS2A_aquatic_screening_plot_table.tsv	2520774	044b020c072966f321ad71d236ff52dcb38a83b59c77992ad479841332ee24c5	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:30:source_data/FigS2_aquatic_single_gene_screen/endpointfix_FigS2B_aquatic_asymmetry_summary.tsv	2675	0e25bedd2e4a45ac379ec716bb841da8291b1d2c8db9e5f3e2e191c1d4806d2e	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:31:source_data/FigS2_aquatic_single_gene_screen/endpointfix_FigS2C_aquatic_enrichment_terms_for_plot.tsv	450291	6056bed74a6cdab279c97278107ca9a2bb572622b1c879dbaf72e231714f4c14	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:32:source_data/Foldwise_ASR_sensitivity/foldwise_ASR_candidate_gene_count_summary.tsv	585	4ec919eca1ff473d5528610b8dd8a1a51b4e6e8d5532f9ed03a26bbad02254b9	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:33:source_data/Foldwise_ASR_sensitivity/foldwise_ASR_fold_level_summary.tsv	214586	13ae3e1c9793aa9c7dd45376bfca59ed6a286d5241420de905866be18e1103d8	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:34:source_data/Foldwise_ASR_sensitivity/foldwise_ASR_selected_predictor_count_summary.tsv	583	761ab5eefa0e81b266467524a314ed27583b28386d8024273b89cc7191695bdc	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./data_manifest/source_data_manifest.tsv:35:source_data/Foldwise_ASR_sensitivity/foldwise_ASR_vs_frozen_ASR_AUC_summary.tsv	1675	7231332bb0d1baf5b79382511b007011305e962d088ef2854039099395902fd5	Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present
./code/05_GBI_construction/build_endpointfix_gbi_matrix.R:56:classify_token <- function(x) {
./code/05_GBI_construction/build_endpointfix_gbi_matrix.R:57:  state_tokens <- c("NA", "NA_fuse", "NA_struct", "NA_topo", "residual_NA", "NaN", "Inf", "-Inf", "")
./code/05_GBI_construction/build_endpointfix_gbi_matrix.R:59:  out[x %in% state_tokens] <- x[x %in% state_tokens]
./code/05_GBI_construction/build_endpointfix_gbi_matrix.R:62:  out[is.na(num) & !(x %in% state_tokens)] <- "unexpected_state"
./code/05_GBI_construction/build_endpointfix_gbi_matrix.R:92:state_matrix <- matrix(classify_token(as.vector(value_chr)), nrow = nrow(value_chr), ncol = ncol(value_chr))
./code/05_GBI_construction/build_endpointfix_gbi_matrix.R:98:non_numeric_tokens <- c("NA", "NA_fuse", "NA_struct", "NA_topo", "residual_NA", "NaN", "Inf", "-Inf", "empty")
./code/05_GBI_construction/build_endpointfix_gbi_matrix.R:99:value_chr[value_chr %in% non_numeric_tokens] <- NA
./code/11_figure_generation/Figure4/prepare_Figure4C_display_tables_final.py:69:    "secretion. Use only as cautious context in Supplementary Table S5 or Discussion."
./code/11_figure_generation/Figure4/prepare_Figure4C_display_tables_final.py:298:            "not evidence for direct altered melatonin secretion."
./env/R_package_version_summary.tsv:399:tokenizers	0.3.0	/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library	NA
./.gitignore:41:secrets*
./.gitignore:47:rejection*
```

Interpretation: no private links, credentials, correspondence, AI chat logs or submission/rejection documents were detected. Hits shown above are false positives caused by biological words containing `secret`, code token-classification variable names, the CRAN package `tokenizers`, or `.gitignore` deny-list patterns.

## Local absolute path check

Command:
```bash
rg -n "/Users/jiaqiwu|/Users/|/Volumes/|/home/|/mnt/" . --glob "!.git/**" --glob "!docs/release_QC_report.md" || true
```

Output:
```text
```

Interpretation: PASS; no local absolute paths remain in the public release tree outside the QC command text itself.

## Large-file check

Command:
```bash
find . -type f -size +100M -not -path "./.git/*" -print
```

Output:
```text
```

Interpretation: PASS; no files larger than 100 MB are included.

## Smoke tests

```text
During startup - Warning messages:
1: Setting LC_CTYPE failed, using "C" 
2: Setting LC_COLLATE failed, using "C" 
3: Setting LC_TIME failed, using "C" 
4: Setting LC_MESSAGES failed, using "C" 
5: Setting LC_MONETARY failed, using "C" 
PASS source data files present
During startup - Warning messages:
1: Setting LC_CTYPE failed, using "C" 
2: Setting LC_COLLATE failed, using "C" 
3: Setting LC_TIME failed, using "C" 
4: Setting LC_MESSAGES failed, using "C" 
5: Setting LC_MONETARY failed, using "C" 
PASS trait table row count = 302
PASS Fig5B public values found
PASS forbidden legacy values absent from public-facing staged assets
```

Interpretation: PASS. Locale warnings from R are environment warnings only.

## Git status before commit

```text
A  .gitignore
A  CITATION.cff
A  LICENSE
A  README.md
A  code/01_alignment_QC/Alignment_length.pl
A  code/01_alignment_QC/CodonAlign.pl
A  code/01_alignment_QC/Concatenate.pl
A  code/01_alignment_QC/CopyFile.pl
A  code/01_alignment_QC/FragmentCount.pl
A  code/01_alignment_QC/FragmentScreen.pl
A  code/01_alignment_QC/MissingSeqCountOnly.pl
A  code/01_alignment_QC/MoveFasta.pl
A  code/01_alignment_QC/MoveFile.pl
A  code/01_alignment_QC/PickUpGenes.pl
A  code/01_alignment_QC/QC_compare/after_QC/A1BG.afterQC.species.txt
A  code/01_alignment_QC/QC_compare/before_QC/A1BG.beforeQC.species.txt
A  code/01_alignment_QC/README.md
A  code/01_alignment_QC/RandomSamplingSites.pl
A  code/01_alignment_QC/RandomSamplingSites_codon.pl
A  code/01_alignment_QC/RemoveEmptySequences.pl
A  code/01_alignment_QC/RemoveFragmentSeq_cds.pl
A  code/01_alignment_QC/RemoveGap_cds.pl
A  code/01_alignment_QC/RemoveGap_cds_length.pl
A  code/01_alignment_QC/RemoveSmallFrag.pl
A  code/01_alignment_QC/SeparateData.pl
A  code/01_alignment_QC/Standard.txt
A  code/01_alignment_QC/Translate.pl
A  code/01_alignment_QC/change_gene_species_names.pl
A  code/01_alignment_QC/hg38.option.sh
A  code/01_alignment_QC/option.sh
A  code/01b_species_missing_gene_QC/6366Genes.txt
A  code/01b_species_missing_gene_QC/Alignment_length.pl
A  code/01b_species_missing_gene_QC/CopyFasta.pl
A  code/01b_species_missing_gene_QC/Data.R
A  code/01b_species_missing_gene_QC/Fasta2phylic.pl
A  code/01b_species_missing_gene_QC/MissingGene40.txt
A  code/01b_species_missing_gene_QC/MissingSeqCountOnly.pl
A  code/01b_species_missing_gene_QC/MoveFasta.pl
A  code/01b_species_missing_gene_QC/PickUpGenes.pl
A  code/01b_species_missing_gene_QC/ScreeningData.R
A  code/01b_species_missing_gene_QC/SeparateData.pl
A  code/01b_species_missing_gene_QC/SeparateData_phy.pl
A  code/01b_species_missing_gene_QC/Sequence_Info.txt
A  code/01b_species_missing_gene_QC/Species_sum.pl
A  code/01b_species_missing_gene_QC/allGenes17434.txt
A  code/01b_species_missing_gene_QC/all_tree.2249.success.txt
A  code/01b_species_missing_gene_QC/failed.623.txt
A  code/01b_species_missing_gene_QC/failed137.txt
A  code/01b_species_missing_gene_QC/failed26.length.txt
A  code/01b_species_missing_gene_QC/failed26.list.txt
A  code/01b_species_missing_gene_QC/gene_qcPass.singleCopy.2275.geneList.txt
A  code/01b_species_missing_gene_QC/hg38_species_name_cleaned.all_species.txt
A  code/01b_species_missing_gene_QC/hg38_species_name_cleaned.length.txt
A  code/01b_species_missing_gene_QC/mammal302.noMissing.geneList.476.txt
A  code/01b_species_missing_gene_QC/mammal302.noMissing_singleCopy.176.geneList.txt
A  code/01b_species_missing_gene_QC/mammal377.MissingCount.txt
A  code/01b_species_missing_gene_QC/mammal377.MissingInfo.txt
A  code/01b_species_missing_gene_QC/mammal377.SpeciesMissingCount.txt
A  code/01b_species_missing_gene_QC/mammal377.gene_qcPass.singleCopy.2310.geneList.txt
A  code/01b_species_missing_gene_QC/mammal377.species.txt
A  code/01b_species_missing_gene_QC/seq_length.txt
A  code/01b_species_missing_gene_QC/species.txt
A  code/01b_species_missing_gene_QC/species301.MissingCount.txt
A  code/01b_species_missing_gene_QC/species301.MissingInfo.txt
A  code/01b_species_missing_gene_QC/species301.SpeciesMissingCount.txt
A  code/01b_species_missing_gene_QC/species301.txt
A  code/01b_species_missing_gene_QC/species302.MissingCount.txt
A  code/01b_species_missing_gene_QC/species302.MissingInfo.txt
A  code/01b_species_missing_gene_QC/species302.SpeciesMissingCount.txt
A  code/01b_species_missing_gene_QC/species_keep302.txt
A  code/02_species_tree_iqtree/_partition.txt
A  code/02_species_tree_iqtree/autoIQTree_nofolder.pl
A  code/02_species_tree_iqtree/loop_run_iqtree.pl
A  code/02_species_tree_iqtree/seq_length.txt
A  code/03_paml_baseml_fixed_topology/AutoBaseML_phy.pl
A  code/03_paml_baseml_fixed_topology/ExtractBaseMLResults2.pl
A  code/03_paml_baseml_fixed_topology/Fasta2phylip.pl
A  code/03_paml_baseml_fixed_topology/SubTree.R
A  code/03_paml_baseml_fixed_topology/_baseml.ctl
A  code/04_SplitAligner/CITATION.cff
A  code/04_SplitAligner/LICENSE
A  code/04_SplitAligner/README.md
A  code/04_SplitAligner/SplitAligner.pl
A  code/04_SplitAligner/scripts/classify_fix_missingness.pl
A  code/04_SplitAligner/scripts/confirm_na_structure.pl
A  code/04_SplitAligner/scripts/extract_na_fuse.pl
A  code/04_SplitAligner/scripts/generate_branch_matrix.pl
A  code/04_SplitAligner/scripts/label_species_tree.pl
A  code/04_SplitAligner/scripts/split_branch_label.pl
A  code/04_SplitAligner/scripts/tree_to_splits.pl
A  code/05_GBI_construction/build_endpointfix_gbi_matrix.R
A  code/06_global_single_gene_screen/build_stage1_ttest_qc.R
A  code/06_global_single_gene_screen/marine_functions.R
A  code/06_global_single_gene_screen/run_ttest_screening.R
A  code/06b_deterministic_ASR_labels/build_stage1_deterministic_ttest_qc.R
A  code/06b_deterministic_ASR_labels/deterministic_asr_precheck.R
A  code/06b_deterministic_ASR_labels/marine_functions.R
A  code/06b_deterministic_ASR_labels/plot_stage1_ancestral_states.R
A  code/06b_deterministic_ASR_labels/run_ttest_screening.R
A  code/06c_positive_count_permutation_control/run_endpointfix_drop_whale_permutation_control.R
A  code/11_figure_generation/Figure1/Clade_static.R
A  code/11_figure_generation/Figure1/Figure1_ggplot.R
A  code/11_figure_generation/Figure1/ggtree_example.R
A  code/11_figure_generation/Figure3/make_Figure3A_endpointfix_volcano.R
A  code/11_figure_generation/Figure3/make_Figure3B_endpointfix_asymmetry.R
A  code/11_figure_generation/Figure3/make_Figure3C_endpointfix_human_phenotype.R
A  code/11_figure_generation/Figure3/make_Figure3_endpointfix.R
A  code/11_figure_generation/Figure4/assemble_Figure4_corrected_old_design.py
A  code/11_figure_generation/Figure4/make_Figure4B_coefficient_architecture_corrected.py
A  code/11_figure_generation/Figure4/make_Figure4C_artwork_draft_from_display_tables.py
A  code/11_figure_generation/Figure4/make_Figure4C_module_partition_corrected.py
A  code/11_figure_generation/Figure4/prepare_Figure4C_display_tables_final.py
A  code/11_figure_generation/Figure4/prepare_Figure4C_display_tables_from_pro_annotations.py
A  code/11_figure_generation/FigureS1/run_aquatic_5d_pca.R
A  code/11_figure_generation/FigureS2/make_FigureS2A_endpointfix_aquatic_volcano.R
A  code/11_figure_generation/FigureS2/make_FigureS2B_endpointfix_aquatic_asymmetry.R
A  code/11_figure_generation/FigureS2/make_FigureS2C_endpointfix_aquatic_human_phenotype.R
A  config/example_paths.yaml
A  data/README.md
A  data_manifest/Dryad_large_data_manifest.tsv
A  data_manifest/external_data_sources.tsv
A  data_manifest/source_data_manifest.tsv
A  docs/data_code_availability_draft.md
A  docs/dryad_manifest.md
A  docs/expected_outputs.md
A  docs/known_limitations.md
A  docs/public_release_cleanup_report.tsv
AM docs/release_QC_report.md
A  docs/release_plan.md
A  docs/reproduce_map.tsv
A  docs/reproduction_guide.md
A  docs/workflow_overview.md
A  env/README_environment_captures.md
A  env/R_package_version_summary.tsv
A  env/R_sessionInfo.txt
A  env/R_sessionInfo_current_local.txt
A  env/python_pip_freeze.txt
A  metadata/branch_label_mapping/README_branch_label_mapping.md
A  metadata/branch_label_mapping/endpointfix_branch_label_crosswalk_new_to_old.tsv
A  metadata/branch_label_mapping/endpointfix_branch_label_crosswalk_old_to_new.tsv
A  metadata/branch_label_mapping/endpointfix_branch_label_crosswalk_summary.tsv
A  metadata/branch_label_mapping/endpointfix_new_branch_split_keys.tsv
A  metadata/branch_label_mapping/endpointfix_old_branch_split_keys.tsv
A  provenance_scripts/README.md
A  provenance_scripts/nonportable_script_manifest.tsv
A  reproduce_map.tsv
A  scripts/README.md
A  session_info.txt
A  source_data/Baseline_slow_gene_layer_overlap/Supplementary_Table_S4_baseline_slow_gene_layer_overlap.tsv
A  source_data/Baseline_slow_gene_layer_overlap/slow_gene_direction_crosstab.tsv
A  source_data/Baseline_slow_gene_layer_overlap/slow_gene_overlap_partitions.tsv
A  source_data/Baseline_slow_gene_layer_overlap/slow_gene_set_counts.tsv
A  source_data/Fig1_trait_framework/Figure1C_score_clade_counts.tsv
A  source_data/Fig1_trait_framework/trait_table.mammal302.active_TY_NK_final_18pt.tsv
A  source_data/Fig2_nested_gLOOCV/Figure2B_nested_ttest_ROC_data.tsv
A  source_data/Fig2_nested_gLOOCV/Figure2C_nested_ttest_OOF_distribution_data.tsv
A  source_data/Fig2_nested_gLOOCV/Figure2_nested_ttest_plot_summary.tsv
A  source_data/Fig3_marine_single_gene_screen/endpointfix_Fig3A_screening_plot_table.tsv
A  source_data/Fig3_marine_single_gene_screen/endpointfix_Fig3B_directional_asymmetry_summary.tsv
A  source_data/Fig3_marine_single_gene_screen/endpointfix_Fig3C_enrichment_terms_for_plot.tsv
A  source_data/Fig4_corrected_full_data_architecture/Figure4A_plot_table.tsv
A  source_data/Fig4_corrected_full_data_architecture/Figure4B_plot_table.tsv
A  source_data/Fig4_corrected_full_data_architecture/Figure4C_TableS5_annotation_table_final.tsv
A  source_data/Fig4_corrected_full_data_architecture/Figure4C_main_display_cell_summary_final.tsv
A  source_data/Fig4_corrected_full_data_architecture/Figure4C_plot_table.tsv
A  source_data/Fig5_sensitivity_permutation_turnover/Figure5A_nested_sensitivity_plotdata.tsv
A  source_data/Fig5_sensitivity_permutation_turnover/Figure5A_nested_sensitivity_table.tsv
A  source_data/Fig5_sensitivity_permutation_turnover/Figure5B_endpointfix_positive_count_permutation_plotdata.tsv
A  source_data/Fig5_sensitivity_permutation_turnover/Figure5B_endpointfix_source_check.tsv
A  source_data/Fig5_sensitivity_permutation_turnover/Figure5C_cross_layer_comparison_supplementary.tsv
A  source_data/Fig5_sensitivity_permutation_turnover/Figure5C_module_null_summary_2col.tsv
A  source_data/Fig5_sensitivity_permutation_turnover/Figure5C_predictor_turnover_metrics_2col.tsv
A  source_data/FigS1_aquaticity_PCA/FigureS1_PCA_plot_tables.tsv
A  source_data/FigS1_aquaticity_PCA/FigureS1_PCA_summary.txt
A  source_data/FigS1_aquaticity_PCA/aquaticity_5d_scores.mammal302.slim_TY_NK_final_18pt.tsv
A  source_data/FigS2_aquatic_single_gene_screen/endpointfix_FigS2A_aquatic_screening_plot_table.tsv
A  source_data/FigS2_aquatic_single_gene_screen/endpointfix_FigS2B_aquatic_asymmetry_summary.tsv
A  source_data/FigS2_aquatic_single_gene_screen/endpointfix_FigS2C_aquatic_enrichment_terms_for_plot.tsv
A  source_data/Foldwise_ASR_sensitivity/foldwise_ASR_candidate_gene_count_summary.tsv
A  source_data/Foldwise_ASR_sensitivity/foldwise_ASR_fold_level_summary.tsv
A  source_data/Foldwise_ASR_sensitivity/foldwise_ASR_selected_predictor_count_summary.tsv
A  source_data/Foldwise_ASR_sensitivity/foldwise_ASR_vs_frozen_ASR_AUC_summary.tsv
A  source_data/README.md
A  tests/smoke_tests/check_expected_counts.R
A  tests/smoke_tests/check_forbidden_old_values.sh
A  tests/smoke_tests/check_source_data_files.R
A  tests/smoke_tests/expected_endpointfix_values.tsv
?? code/07_nested_ttest_baseline_gLOOCV/
?? code/08_Figure4_Figure5_alignment/
?? code/09_corrected_full_data_LASSO_architecture/
?? code/10_foldwise_ASR_sensitivity/
?? code/11_figure_generation/Figure2/
?? code/11_figure_generation/Figure4/prepare_inputs_and_make_Figure4A_corrected.R
?? code/11_figure_generation/Figure5/
```
